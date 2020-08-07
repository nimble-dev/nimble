## Used when syntax xi[1:N] ~ dCRP(conc) is used in BUGS.

## FIXME: this is temporary function until we bring this into the full model API
getParentNodes <- function(nodes, model, returnType = 'names', stochOnly = FALSE) {
    ## adapted from BUGS_modelDef creation of edgesFrom2To
    getParentNodesCore <- function(nodes, model, returnType = 'names', stochOnly = FALSE) {
        nodeIDs <- model$expandNodeNames(nodes, returnType = "ids")
        fromIDs <- sort(unique(unlist(edgesTo2From[nodeIDs])))
        fromNodes <- maps$graphID_2_nodeName[fromIDs]
        if(!length(fromNodes))
            return(character(0))
        fromNodesDet <- fromNodes[model$modelDef$maps$types[fromIDs] == 'determ']
        ## Recurse through parents of deterministic nodes.
        fromNodes <- c(if(stochOnly) fromNodes[model$modelDef$maps$types[fromIDs] == 'stoch'] else fromNodes, 
                       if(length(fromNodesDet)) getParentNodesCore(fromNodesDet, model, returnType, stochOnly) else character(0))
        fromNodes
    }
    
    maps <- model$modelDef$maps
    maxNodeID <- length(maps$vertexID_2_nodeID) ## should be same as length(maps$nodeNames)
    ## Only determine edgesTo2From once and then obtain in getParentNodesCore via scoping.
    edgesLevels <- if(maxNodeID > 0) 1:maxNodeID else numeric(0)
    fedgesTo <- factor(maps$edgesTo, levels = edgesLevels) ## setting levels ensures blanks inserted into the splits correctly
    edgesTo2From <- split(maps$edgesFrom, fedgesTo)

    getParentNodesCore(nodes, model, returnType, stochOnly)
}

getNsave <- nimbleFunction(
  setup = function(mvSaved){
    niter <- 0
  },
  run = function(){
    niter <<- getsize(mvSaved) # number of iterations in the MCMC
  },
  methods = list( reset = function () {} )
)

getParentNodes <- function(nodes, model, returnType = 'names', stochOnly = FALSE) {
  ## adapted from BUGS_modelDef creation of edgesFrom2To
  getParentNodesCore <- function(nodes, model, returnType = 'names', stochOnly = FALSE) {
    nodeIDs <- model$expandNodeNames(nodes, returnType = "ids")
    fromIDs <- sort(unique(unlist(edgesTo2From[nodeIDs])))
    fromNodes <- maps$graphID_2_nodeName[fromIDs]
    if(!length(fromNodes))
      return(character(0))
    fromNodesDet <- fromNodes[model$modelDef$maps$types[fromIDs] == 'determ']
    ## Recurse through parents of deterministic nodes.
    fromNodes <- c(if(stochOnly) fromNodes[model$modelDef$maps$types[fromIDs] == 'stoch'] else fromNodes, 
                   if(length(fromNodesDet)) getParentNodesCore(fromNodesDet, model, returnType, stochOnly) else character(0))
    fromNodes
  }
  
  maps <- model$modelDef$maps
  maxNodeID <- length(maps$vertexID_2_nodeID) ## should be same as length(maps$nodeNames)
  ## Only determine edgesTo2From once and then obtain in getParentNodesCore via scoping.
  edgesLevels <- if(maxNodeID > 0) 1:maxNodeID else numeric(0)
  fedgesTo <- factor(maps$edgesTo, levels = edgesLevels) ## setting levels ensures blanks inserted into the splits correctly
  edgesTo2From <- split(maps$edgesFrom, fedgesTo)
  
  getParentNodesCore(nodes, model, returnType, stochOnly)
}

##-----------------------------------------
##  Wrapper function for sampleDPmeasure
##-----------------------------------------

#' Get posterior samples for a Dirichlet process measure
#'
#' This function obtains posterior samples from a Dirichlet process distributed random measure of a model specified using the \code{dCRP} distribution.
#'
#' @param MCMC an MCMC class object, either compiled or uncompiled.
#' @param epsilon  used for determining the truncation level of the representation of the random measure.
#' @param setSeed Logical or numeric argument. If a single numeric value is provided, R's random number seed will be set to this value. In the case of a logical value, if \code{TRUE}, then R's random number seed will be set to \code{1}. Note that specifying the argument \code{setSeed = 0} does not prevent setting the RNG seed, but rather sets the random number generation seed to \code{0}.  Default value is \code{FALSE}.
#' 
#' @author Claudia Wehrhahn and Christopher Paciorek
#' 
#' @export
#' @details
#' This function provides samples from a random measure having a Dirichlet process prior. Realizations are almost surely discrete and represented by a (finite) stick-breaking representation (Sethuraman, 1994), whose atoms (or point masses) are independent and identically distributed. This sampler can only be used with models containing a \code{dCRP} distribution. 
#'
#' The \code{MCMC} argument is an object of class MCMC provided by \code{buildMCMC}, or its compiled version. The MCMC should already have been run, as \code{getSamplesDPmeasure} uses the posterior samples to generate samples of the random measure. Note that the monitors associated with that MCMC must include the cluster membership variable (which has the \code{dCRP} distribution), the cluster parameter variables, all variables directly determining the \code{dCRP} concentration parameter, and any stochastic parent variables of the cluster parameter variables. See \code{help(configureMCMC)} or \code{help(addMonitors)} for information on specifying monitors for an MCMC.
#' 
#' The \code{epsilon} argument is optional and used to determine the truncation level of the random measure. \code{epsilon} is the tail probability of the random measure, which together with posterior samples of the concentration parameter, determines the truncation level. The default value is 1e-4.
#'  
#' The output is a list of matrices. Each matrix represents a sample from the random measure. In order to reduce the output's dimensionality, the weigths of identical atoms are added up. The stick-breaking weights are named \code{weights} and the atoms are named based on the cluster variables in the model.
#' 
#' For more details about sampling the random measure and determining its truncation level, see Section 3 in Gelfand, A.E. and Kottas, A. 2002.
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
#' }
getSamplesDPmeasure <- function(MCMC, epsilon = 1e-4, setSeed = FALSE) {
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
  
  rsampler <- sampleDPmeasure(model, mvSamples, epsilon) # 

  niter <- getsize(MCMC$mvSamples)
  samplesMeasure <- list()
  
  if(is.numeric(setSeed)) {
    set.seed(setSeed[1])
    if(length(setSeed) > 1) {
      nimCat('getSamplesDPmeasure: setSeed argument has length > 1 and only the first element will be used') 
    }
  } else if(setSeed) set.seed(1)
  
  if(compiled) {
    csampler <- compileNimble(rsampler, project = model)
    for(i in 1:niter) {
      samplesMeasure[[i]] <- csampler$run(i)
    }
  } else {
    for(i in 1:niter) {
      samplesMeasure[[i]] <- rsampler$run(i)
    }
  }
  
  
  dcrpVar <- rsampler$dcrpVar
  clusterVarInfo <- findClusterNodes(model, dcrpVar) 
  p <- rsampler$p
  namesVars <- getSamplesDPmeasureNames(clusterVarInfo, model, 1, p)
  for(i in 1:niter) {
    colnames(samplesMeasure[[i]]) <- c("weights", namesVars)
  }
  
  return(samplesMeasure)
}


sampleDPmeasure <- nimbleFunction(
  name = 'sampleDPmeasure',

  setup = function(model, mvSaved, epsilon) { # 
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
    
    ## Check that cluster parameters are IID (across clusters, non-IID within clusters is ok), as required for random measure G
    isIID <- TRUE
    for(i in seq_along(clusterVarInfo$clusterNodes)) {
      clusterNodes <- clusterVarInfo$clusterNodes[[i]]  # e.g., 'thetatilde[1]',...,
      clusterIDs <- clusterVarInfo$clusterIDs[[i]]
      splitNodes <- split(clusterNodes, clusterIDs)
      valueExprs <- lapply(splitNodes, function(x) {
        out <- sapply(x, model$getValueExpr)
        names(out) <- NULL
        out
      })
      if(length(unique(valueExprs)) != 1) 
        isIID <- FALSE
    }
    
    if(!isIID && length(tildeVars) == 2 && checkNormalInvGammaConjugacy(model, clusterVarInfo, length(dcrpElements), 'dinvgamma'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && checkNormalInvGammaConjugacy(model, clusterVarInfo, length(dcrpElements), 'dgamma'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && checkNormalInvWishartConjugacy(model, clusterVarInfo, length(dcrpElements), 'dinvwish'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && checkNormalInvWishartConjugacy(model, clusterVarInfo, length(dcrpElements), 'dwish'))
      isIID <- TRUE
    ## Tricky as MCMC might not be using conjugacy, but presumably ok to proceed regardless of how
    ## MCMC was done, since conjugacy existing would guarantee IID.
    if(!isIID) stop('sampleDPmeasure: cluster parameters have to be independent and identically distributed. \n')
    
    ## Check that necessary variables are being monitored.
    
    ## Check that cluster variables are monitored.
    counts <- tildeVars %in% mvSavedVars
    if( sum(counts) != length(tildeVars) ) 
      stop('sampleDPmeasure: The node(s) representing the cluster variables must be monitored in the MCMC (and therefore stored in the modelValues object).\n')  
    
    parentNodesTildeVars <- getParentNodes(tildeVars, model, returnType = 'names', stochOnly = TRUE)
    if(length(parentNodesTildeVars)) {
      parentNodesTildeVarsDeps <- model$getDependencies(parentNodesTildeVars, self = FALSE)
    } else parentNodesTildeVarsDeps <- NULL
    ## make sure tilde nodes are included (e.g., if a tilde node has no stoch parents) so they get simulated
    parentNodesTildeVarsDeps <- model$topologicallySortNodes(c(parentNodesTildeVarsDeps, unlist(clusterVarInfo$clusterNodes)))
    
    if(!all(model$getVarNames(nodes = parentNodesTildeVars) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the cluster variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')

    ## Check that parent nodes of cluster IDs are monitored.  
    parentNodesXi <- getParentNodes(dcrpNode, model, returnType = 'names', stochOnly = TRUE)
    
    if(!all(model$getVarNames(nodes = parentNodesXi) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the membership variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    
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
    N <- length(model$expandNodeNames(dcrpNode, returnScalarComponents = TRUE))
    
    p <- length(tildeVars)
    lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE))
    dimTildeVarsNim <- numeric(p+1) # nimble dimension (0 is scalar, 1 is 2D array, 2 is 3D array) (dimTildeVarsNim=dimTildeNim)
    dimTildeVars <- numeric(p+1) # dimension to be used in run code (dimTildeVars=dimTilde)
    for(i in 1:p) {
      dimTildeVarsNim[i] <- model$getDimension(clusterVarInfo$clusterNodes[[i]][1])
      dimTildeVars[i] <- lengthData^(dimTildeVarsNim[i]) 
    }
    nTildeVarsPerCluster <-  clusterVarInfo$numNodesPerCluster
    nTilde <- numeric(p+1)
    nTilde[1:p] <- clusterVarInfo$nTilde / nTildeVarsPerCluster
    if(any(nTilde[1:p] != nTilde[1])){
      stop('sampleDPmeasure: All cluster parameters must have the same number of parameters.\n')
    }
    
    tildeVarsCols <- c(dimTildeVars[1:p]*nTildeVarsPerCluster, 0)
    tildeVarsColsSum <- c(0, cumsum(tildeVarsCols))
    mvIndexes <- matrix(0, nrow=nTilde[1], ncol=(sum(dimTildeVars[1:p]*nTildeVarsPerCluster))) 
    for(j in 1:p) {
      tildeNodesModel <- model$expandNodeNames(clusterVarInfo$clusterVars[j], returnScalarComponents=TRUE) # tilde nodes j in model
      allIndexes <- 1:length(tildeNodesModel)
      for(l in 1:nTilde[1]) {
        clusterID <- l
        tildeNodesPerClusterID <- model$expandNodeNames(clusterVarInfo$clusterNodes[[j]][clusterVarInfo$clusterIDs[[j]] == clusterID], returnScalarComponents=TRUE) # tilde nodes in cluster with id 1
        aux <- match(tildeNodesModel, tildeNodesPerClusterID, nomatch = 0) 
        mvIndexes[l,(tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1] ] <- which(aux != 0)
      }
    }
    
    ## control list extraction
    ## The error of approximation G is given by (conc / (conc +1))^{truncG-1}. 
    ## we are going to define an error of approximation and based on the posterior values of the conc parameter define the truncation level of G
    ## the error is between errors that are considered very very small in the folowing papers
    ## Ishwaran, H., & James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. Journal of the American Statistical Association, 96(453), 161-173.
    ## Ishwaran, H., & Zarepour, M. (2000). Markov chain Monte Carlo in approximate Dirichlet and beta two-parameter process hierarchical models. Biometrika, 87(2), 371-390.
    # epsilon <- 1e-4
    setupOutputs(dcrpVar)
    
  },
  run = function(m = integer()) {
    returnType(double(2))
    
    ## value of concentraton parameter
    if( fixedConc ) { 
      concSamples <- model$getParam(dcrpNode, 'conc')
    } else { 
      nimCopy(from = mvSaved, to = model, nodes = parentNodesXi, row = m) 
      model$calculate(parentNodesXiDeps)
      concSamples <- model$getParam(dcrpNode, 'conc')
    }
    
    ## truncation level
    concAux <- concSamples + N
    truncG <- log(epsilon) / log(concAux / (concAux+1)) 
    truncG <- ceiling(truncG) 
    
    ## getting the unique cluster variables that where sampled in the mcmc and the sampling probabilities of the polya urn of the unique cluster variables
    probs <- nimNumeric(N) # polya urn probabilities
    uniqueValues <- matrix(0, nrow = N, ncol = tildeVarsColsSum[p+1])  # unique cluster variables
    xiiter <- mvSaved[dcrpVar, m]
    range <- min(xiiter):max(xiiter) 
    index <- 1
    for(i in seq_along(range)){   
      cond <- sum(xiiter == range[i])
      if(cond > 0){
        probs[index] <- cond
        ## slight workaround because can't compile mvSaved[tildeVars[j], m]  
        nimCopy(mvSaved, model, tildeVars, row = m)
        for(j in 1:p){
          jcols <- (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]
          uniqueValues[index, jcols] <- values(model, tildeVars[j])[mvIndexes[range[i], jcols]]
        }
        index <- index+1
      }
    }
    probs[index] <- concSamples # last probability of the polya urn (the sample concetration parameter)
    newValueIndex <- index 
    
    ## copy tilde parents into model for use in simulation below when simulate atoms of G_0  
    nimCopy(mvSaved, model, parentNodesTildeVars, row = m)
    
    
    ## sampling random measure:
    # this is the reduced version in the sense that weights of identical atoms are added up
    # weights and atom are sampled in the same loop
    
    ## sampling stick-breaking weights: as many as truncG
        
    ## reduced samples from random measure: 
    weights <- nimNumeric(truncG) # weights of random measure
    weightsTmp <- nimNumeric(truncG)
    atoms <- matrix(0, ncol = tildeVarsColsSum[p+1], nrow = truncG) # atoms of random measure
    indexesG <- nimNumeric(truncG) # indicates if an existing or a new atom is sampled from the polya urn. New tom are indicated by newValueIndex 
    indexG0 <- newValueIndex # used for new atoms
    
    #vaux <- rbeta(1, 1, concAux)
    #v1prod <- 1
    #weightsTmp[1] <- vaux  
    #for(l1 in 2:truncG) {
    #  v1prod <- v1prod * (1-vaux)
    #  vaux <- rbeta(1, 1, concAux)
    #  weightsTmp[l1] <- vaux * v1prod 
    #}
    #weightsTmp[1:truncG] <- weightsTmp[1:truncG] / (1 - v1prod * (1-vaux)) # normalizing weigths
    
    vaux <- rbeta(1, 1, concAux)
    v1prod <- 1
    for(l1 in 1:truncG) {
      # sampling a weight
      if(l1 == 1) {
        weightsTmp[l1] <- vaux 
      } else {
        v1prod <- v1prod * (1-vaux)
        vaux <- rbeta(1, 1, concAux)
        weightsTmp[l1] <- vaux * v1prod 
      }

      # sampling an atom
      indexesG[l1] <- rcat(prob = probs[1:newValueIndex])
      if(indexesG[l1] < newValueIndex) { # an existing atom was sampled and the corresponding weights are added
        weights[indexesG[l1]] <- weights[indexesG[l1]] + weightsTmp[l1] 
        sumCol <- 0
        for(j in 1:p){
          jcols <- (sumCol + 1):(sumCol + tildeVarsCols[j]) 
          atoms[indexesG[l1], jcols] <-  uniqueValues[indexesG[l1], (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]]  
          sumCol <- sumCol + tildeVarsCols[j] 
        }
      } else { # a new atom is sampled (from G0) with its corresponding new weight
        weights[indexG0] <-  weightsTmp[l1] 
        model$simulate(parentNodesTildeVarsDeps)
        sumCol <- 0
        for(j in 1:p){
          jcols <- (sumCol + 1):(sumCol + tildeVarsCols[j])
          atoms[indexG0, jcols] <- values(model, tildeVars[j])[mvIndexes[1,  (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]]] # <<-
          sumCol <- sumCol + tildeVarsCols[j] 
        }
        indexG0 <- indexG0 + 1
      }
    }
    weights[1:truncG] <- weights[1:truncG] / (1 - v1prod * (1-vaux)) # normalizing weigths

    ## check that all unique tilde variables were actually sampled. If not we need to rearrange when creating final output
    missingIndex <- nimNumeric(newValueIndex-1)
    uniqueIndex <- 1:(newValueIndex - 1)
    ii <- 1
    for(i in seq_along(uniqueIndex)) {
      if( !any(indexesG == uniqueIndex[i]) ) {
        missingIndex[ii] <- uniqueIndex[i]
        ii <- ii + 1
      }
    }
    
    # final sample of the random measure saved in 'output'
    nrowG <- indexG0 - 1 - (ii - 1)
    output <- nimNumeric(nrowG*(tildeVarsColsSum[p+1]+1))
    ii <- 1
    imissing <- 1 
    for(i in 1:(indexG0 - 1)) { # saving weights
      if(i != missingIndex[imissing]) {
        output[ii] <- weights[i]
        ii <- ii + 1  
      } else {
        imissing <- imissing + 1
      }
    }
    for(j in 1:tildeVarsColsSum[p+1]) { # saving atoms
      imissing <- 1
      for(i in 1:(indexG0 - 1)) {
        if(i != missingIndex[imissing]) {
          output[ii] <- atoms[i, j]
          ii <- ii + 1  
        } else {
          imissing <- imissing + 1
        }
      }
    }
    
    outputG <- matrix(nimNumeric(value = output, length=(nrowG*(tildeVarsColsSum[p+1]+1))), ncol=(tildeVarsColsSum[p+1]+1), nrow=nrowG)
    
    return(outputG)
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
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {  ## calculate prior predictive for new cluster for conjugate cases
      returnType(double())
    },
    sample = function(i = integer(), j = integer()) {}  ## sample parameter values for new cluster
  )
)


CRP_nonconjugate <- nimbleFunction(
  name = "CRP_nonconjugate",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    len <- length(model$expandNodeNames(marginalizedNodes[1:M], returnScalarComponents = TRUE))
    saved <- nimNumeric(len+1) # to avoid scalar if M=1
    savedIdx <- 1  
  },
  methods = list(
    storeParams = function() {},  ## nothing needed for non-conjugate
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        out <- out + model$getLogProb(dataNodes[(i-1)*J+j1])
      }
      return(out)
    },
    sample = function(i = integer(), j = integer() ) {
      if(j == 0) {   ## reset to stored values (for case of new cluster not opened)
        values(model, marginalizedNodes[ ((savedIdx-1)*M+1):(savedIdx*M) ]) <<- saved[1:len]
      } else {  ## sample from prior
        savedIdx <<- j  
        saved[1:len] <<- values(model, marginalizedNodes[ ((j-1)*M+1):(j*M) ])
        model$simulate(marginalizedNodes[ ((j-1)*M+1):(j*M) ])
      }
    }
  )
)


## All conjugate samplers assume J==M. This is checked in sampler_CRP.

CRP_conjugate_dnorm_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dnorm_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) { 
    priorMean <- nimNumeric(J+1)
    priorVar <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorMean[j1] <<- model$getParam(marginalizedNodes[j1], 'mean')
        priorVar[j1] <<- model$getParam(marginalizedNodes[j1], 'var')
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        dataVar <- model$getParam(dataNodes[(i-1)*J+j1], 'var')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + dnorm(y, priorMean[j1], sqrt(priorVar[j1] + dataVar), log=TRUE) 
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataVar <- model$getParam(dataNodes[(i-1)*J+j1], 'var')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        postVar <- 1 / (1 / dataVar + 1 / priorVar[j1])
        postMean <- postVar * (y / dataVar + priorMean[j1] / priorVar[j1])
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rnorm(1, postMean, sqrt(postVar)))
      }
    }
  )
)

CRP_conjugate_dmnorm_dmnorm <- nimbleFunction(
  name = "CRP_conjugate_dmnorm_dmnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    d <- length(model[[marginalizedNodes[1]]])
    priorMean <- matrix(0, ncol=d, nrow=J)
      priorCov <- array(0, c(d, d, J))
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorMean[j1, ] <<- model$getParam(marginalizedNodes[j1], 'mean')
        priorCov[, , j1] <<- model$getParam(marginalizedNodes[j1], 'cov') 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        dataCov <- model$getParam(dataNodes[(i-1)*J+j1], 'cov')
        y <- values(model, dataNodes[(i-1)*J+j1])  
        out <- out + dmnorm_chol(y, priorMean[j1, ], chol(priorCov[, , j1] + dataCov), prec_param = FALSE, log=TRUE)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataCov <- model$getParam(dataNodes[(i-1)*J+j1], 'cov')
        y <- values(model, dataNodes[(i-1)*J+j1])
        
        dataPrec <- inverse(dataCov)
        priorPrec <- inverse(priorCov[, , j1])
        postPrecChol <- chol(dataPrec + priorPrec)
        postMean <- backsolve(postPrecChol, forwardsolve(t(postPrecChol), (dataPrec %*% y + priorPrec %*% priorMean[j1, ])[,1]))
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- rmnorm_chol(1, postMean, postPrecChol, prec_param = TRUE) 
      }
    }
  )
)


CRP_conjugate_dnorm_dnorm_nonidentity <- nimbleFunction(
  name = "CRP_conjugate_dnorm_dnorm_nonidentity",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, intermNodes, nInterm, J, M) {
    priorMean <- nimNumeric(J+1)
    priorVar <- nimNumeric(J+1)
    offset <- nimNumeric(J+1)
    coeff <- nimNumeric(J+1)
    currentValue <- nimNumeric(J+1)
    currentInterm <- nimNumeric(nInterm+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
          priorMean[j1] <<- model$getParam(marginalizedNodes[j1], 'mean')
          priorVar[j1] <<- model$getParam(marginalizedNodes[j1], 'var')
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {
        ## In mean of observation, determine a,b in 'a + b*mu[xi[i]]'.
        currentValue <<- values(model, marginalizedNodes[((j-1)*J+1):(j*J)])
        currentInterm <<- values(model, intermNodes[((i-1)*nInterm+1):(i*nInterm)])
        
        values(model, marginalizedNodes[((j-1)*J+1):(j*J)]) <<- rep(0, J)
        model$calculate(intermNodes[((i-1)*nInterm+1):(i*nInterm)])  
        for(j1 in 1:J)
            offset[j1] <<- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        values(model, marginalizedNodes[((j-1)*J+1):(j*J)]) <<- rep(1, J)
        model$calculate(intermNodes[((i-1)*nInterm+1):(i*nInterm)])  
        for(j1 in 1:J)
            coeff[j1] <<- model$getParam(dataNodes[(i-1)*J+j1], 'mean') - offset[j1]
        
        values(model, marginalizedNodes[((j-1)*J+1):(j*J)]) <<- currentValue
        values(model, intermNodes[((i-1)*nInterm+1):(i*nInterm)]) <<- currentInterm
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
          dataVar <- model$getParam(dataNodes[(i-1)*J+j1], 'var')
          y <- values(model, dataNodes[(i-1)*J+j1])[1]
          out <- out + dnorm(y, offset[j1] + coeff[j1]*priorMean[j1], sqrt(coeff[j1]^2 * priorVar[j1] + dataVar), log=TRUE)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
          dataVar <- model$getParam(dataNodes[(i-1)*J+j1], 'var')
          y <- values(model, dataNodes[(i-1)*J+j1])[1]
          postVar <- 1 / (coeff[j1]^2 / dataVar + 1 / priorVar[j1])
          postMean <- postVar * (coeff[j1]*(y-offset[j1]) / dataVar + priorMean[j1] / priorVar[j1])
          values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rnorm(1, postMean, sqrt(postVar)))
      }
    }
  )
)


CRP_conjugate_dnorm_invgamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dnorm_invgamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes1, marginalizedNodes2, dataNodes, J, M) {
    priorMean <- nimNumeric(J+1)
    kappa <- nimNumeric(J+1)
    priorShape <- nimNumeric(J+1)
    priorScale <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorMean[j1] <<- model$getParam(marginalizedNodes1[j1], 'mean')
        kappa[j1] <<- values(model, marginalizedNodes2[j1])[1]/model$getParam(marginalizedNodes1[j1], 'var') # construct kappa as sigma2/(sigma2/kappa)
        priorShape[j1] <<- model$getParam(marginalizedNodes2[j1], 'shape')
        priorScale[j1] <<- model$getParam(marginalizedNodes2[j1], 'scale') 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        c1 <- priorShape[j1] * log(priorScale[j1]) + lgamma(priorShape[j1] + 1/2) + 0.5*log(kappa[j1]) -
          lgamma(priorShape[j1]) - 0.5*log(2) - 0.5*log(pi) - 0.5*log(1 + kappa[j1])
        c2 <- - (priorShape[j1]  + 1/2) * log( (priorScale[j1] + kappa[j1] * (y - priorMean[j1])^2 / (2*(1+kappa[j1])) ) )
        out <- out + c1 + c2
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes2[(j-1)*J+j1]) <<- c(rinvgamma(1, shape = priorShape[j1] + 1/2,
                                                             scale = priorScale[j1] + kappa[j1] * (y - priorMean[j1])^2 / (2*(1+kappa[j1])) ))
        values(model, marginalizedNodes1[(j-1)*J+j1]) <<- c(rnorm(1, mean = (kappa[j1] * priorMean[j1] + y)/(1 + kappa[j1]), 
                                                         sd = sqrt(values(model, marginalizedNodes2[(j-1)*J+j1])[1] / (1+kappa[j1]))) )
      }
    }
  )
)


CRP_conjugate_dinvwish_dmnorm <- nimbleFunction(
  name = "CRP_conjugate_dinvwish_dmnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    d <- length(model[[dataNodes[1]]])
    df0 <- nimNumeric(J+1)
    priorScale <- array(0, c(d, d, J)) # matrix(0, ncol=d, nrow=d)
    c1 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        df0[j1] <<- model$getParam(marginalizedNodes[j1], 'df') 
        priorScale[, , j1] <<- model$getParam(marginalizedNodes[j1], 'S') 
        c1[j1] <<- - 0.5*d*log(2*pi) + 0.5*df0[j1]*logdet(priorScale[, , j1]) + 0.5*d*log(2) +
          sum(lgamma((df0[j1]+2-1:d)/2) - lgamma((df0[j1]+1-1:d)/2)) 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        c2 <- -0.5*(df0[j1]+1) * logdet(priorScale[, , j1] + (y-dataMean)%*%t(y-dataMean) )
        out <- out + c1[j1] + c2
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- rinvwish_chol(1, chol(priorScale[, , j1] + (y-dataMean)%*%t(y-dataMean)),
                                                                       df = (df0[j1]+1), scale_param=TRUE )
      }
    }
  )
)


CRP_conjugate_dwish_dmnorm <- nimbleFunction(
  name = "CRP_conjugate_dwish_dmnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    d <- length(model[[dataNodes[1]]])
    df0 <- nimNumeric(J+1)
    priorRate <- array(0, c(d, d, J)) # matrix(0, ncol=d, nrow=d)
    c1 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        df0[j1] <<- model$getParam(marginalizedNodes[j1], 'df') 
        priorRate[, , j1] <<- model$getParam(marginalizedNodes[j1], 'R') 
        c1[j1] <<- - 0.5*d*log(pi) + 0.5*df0[j1]*logdet(priorRate[, , j1]) +
          sum(lgamma((df0[j1]+2-1:d)/2) - lgamma((df0[j1]+1-1:d)/2)) 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        c2 <- -0.5*(df0[j1]+1) * logdet(priorRate[, , j1] + (y-dataMean)%*%t(y-dataMean) )
        out <- out + c1[j1] + c2
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- rwish_chol(1, chol(priorRate[, , j1] + (y-dataMean)%*%t(y-dataMean)),
                                                                       df = (df0[j1]+1), scale_param=FALSE )
      }
    }
  )
)


CRP_conjugate_dmnorm_invwish_dmnorm <- nimbleFunction(
  name = "CRP_conjugate_dmnorm_invwish_dmnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes1, marginalizedNodes2, dataNodes, J, M) {
    d <- length(model[[marginalizedNodes1[1]]])
    priorMean <- matrix(0, ncol=d, nrow=J)  #nimNumeric(d)
    kappa <- nimNumeric(J+1)
    df0 <- nimNumeric(J+1)
    priorScale <- array(0, c(d, d, J)) # matrix(0, ncol=d, nrow=d)
    c1 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorMean[j1, ] <<- model$getParam(marginalizedNodes1[j1], 'mean') 
        kappa[j1] <<- values(model, marginalizedNodes2[j1])[1] / model$getParam(marginalizedNodes1[j1], 'cov')[1,1]
        df0[j1] <<- model$getParam(marginalizedNodes2[j1], 'df') 
        priorScale[, , j1] <<- model$getParam(marginalizedNodes2[j1], 'S') 
        c1[j1] <<- d*(log(kappa[j1]) - log(1+kappa[j1]) - log(pi))/2 + df0[j1]*logdet(priorScale[, , j1])/2 +
          sum(lgamma((df0[j1]+2-1:d)/2) - lgamma((df0[j1]+1-1:d)/2)) 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        c2 <- -(df0[j1]+1) * logdet(priorScale[, , j1] + (kappa[j1]/(kappa[j1]+1)) * (y-priorMean[j1, ])%*%t(y-priorMean[j1, ]) ) / 2
        out <- out + c1[j1] + c2
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        tmp <- rinvwish_chol(1, chol(priorScale[, , j1] + (kappa[j1]/(kappa[j1]+1)) * (y-priorMean[j1, ])%*%t(y-priorMean[j1, ])),
                             df = (df0[j1]+1), scale_param=TRUE )
        values(model, marginalizedNodes2[(j-1)*J+j1]) <<- c(tmp)
        values(model, marginalizedNodes1[(j-1)*J+j1]) <<- c(rmnorm_chol(1, mean = (kappa[j1] * priorMean[j1, ] + y)/(1 + kappa[j1]), 
                                                               chol( tmp / (1+kappa[j1]) ),
                                                               prec_param = FALSE))
      }
    }
  )
)


## conjugate dnorm_gamma_dnorm
CRP_conjugate_dnorm_gamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dnorm_gamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes1, marginalizedNodes2, dataNodes, J, M) {
    priorMean <- nimNumeric(J+1)
    kappa <- nimNumeric(J+1)
    priorShape <- nimNumeric(J+1)
    priorRate <- nimNumeric(J+1)
    c1 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorMean[j1] <<- model$getParam(marginalizedNodes1[j1], 'mean')
        kappa[j1] <<- model$getParam(marginalizedNodes1[j1], 'tau') / values(model, marginalizedNodes2[j1])[1] # construct kappa
        priorShape[j1] <<- model$getParam(marginalizedNodes2[j1], 'shape')
        priorRate[j1] <<- model$getParam(marginalizedNodes2[j1], 'rate') 
        c1[j1] <<- priorShape[j1] * log(priorRate[j1]) + lgamma(priorShape[j1] + 1/2) + 0.5*log(kappa[j1]) -
          lgamma(priorShape[j1]) - 0.5*log(2) - 0.5*log(pi) - 0.5*log(1 + kappa[j1])
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        c2 <- - (priorShape[j1]  + 1/2) * log( priorRate[j1] + kappa[j1] * (y - priorMean[j1])^2 / (2*(1+kappa[j1]))  )
        out <- out + c1[j1] + c2
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes2[(j-1)*J+j1]) <<- c(rgamma(1, shape = priorShape[j1] + 1/2,
                                                                   rate = priorRate[j1] + kappa[j1] * (y - priorMean[j1])^2 / (2*(1+kappa[j1])) ))
        values(model, marginalizedNodes1[(j-1)*J+j1]) <<- c(rnorm(1, mean = (kappa[j1] * priorMean[j1] + y)/(1 + kappa[j1]), 
                                                                  sd = sqrt(1/(values(model, marginalizedNodes2[(j-1)*J+j1])[1] * (1+kappa[j1]))) ))
      }
    }
  )
)


CRP_conjugate_dmnorm_wish_dmnorm <- nimbleFunction(
  name = "CRP_conjugate_dmnorm_wish_dmnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes1, marginalizedNodes2, dataNodes, J, M) {
    d <- length(model[[marginalizedNodes1[1]]])
    priorMean <- matrix(0, ncol=d, nrow=J)  
    kappa <- nimNumeric(J+1)
    df0 <- nimNumeric(J+1)
    priorRate <- array(0, c(d, d, J)) 
    c1 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorMean[j1, ] <<- model$getParam(marginalizedNodes1[j1], 'mean') 
        kappa[j1] <<- model$getParam(marginalizedNodes1[j1], 'prec')[1,1] / values(model, marginalizedNodes2[j1])[1]
        df0[j1] <<- model$getParam(marginalizedNodes2[j1], 'df') 
        priorRate[, , j1] <<- model$getParam(marginalizedNodes2[j1], 'R') 
        c1[j1] <<- d*(log(kappa[j1]) - log(1+kappa[j1]) - log(pi))/2 + df0[j1]*logdet(priorRate[, , j1])/2 +
          sum(lgamma((df0[j1]+2-1:d)/2) - lgamma((df0[j1]+1-1:d)/2)) 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        c2 <- -(df0[j1]+1) * logdet(priorRate[, , j1] + (kappa[j1]/(kappa[j1]+1)) * (y-priorMean[j1, ])%*%t(y-priorMean[j1, ]) ) / 2
        out <- out + c1[j1] + c2
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        tmp <- rwish_chol(1, chol(priorRate[, , j1] + (kappa[j1]/(kappa[j1]+1)) * (y-priorMean[j1, ])%*%t(y-priorMean[j1, ])),
                             df = (df0[j1]+1), scale_param=FALSE )
        values(model, marginalizedNodes2[(j-1)*J+j1]) <<- c(tmp)
        values(model, marginalizedNodes1[(j-1)*J+j1]) <<- c(rmnorm_chol(1, mean = (kappa[j1] * priorMean[j1, ] + y)/(1 + kappa[j1]), 
                                                                        chol( tmp * (1+kappa[j1]) ),
                                                                        prec_param = TRUE))
      }
    }
  )
)


CRP_conjugate_dinvgamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dinvgamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) { 
    priorShape <- nimNumeric(J+1)
    priorScale <- nimNumeric(J+1)
    c1 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape[j1] <<- model$getParam(marginalizedNodes[j1], 'shape')
        priorScale[j1] <<- model$getParam(marginalizedNodes[j1], 'scale')
        c1[j1] <<- - 0.5 * log(2) - 0.5 * log(pi) + priorShape[j1] * log(priorScale[j1]) - 
          lgamma(priorShape[j1]) + lgamma(priorShape[j1] + 0.5)
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        c2 <- -(priorShape[j1] + 0.5) * log(priorScale[j1] + 0.5 * (y - dataMean)^2)
        out <- out + c1[j1] + c2 
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rinvgamma(1, shape = priorShape[j1]  + 0.5,
                                                                     scale = priorScale[j1] + 0.5 * (y - dataMean)^2))
      }
    }
  )
)






CRP_conjugate_dgamma_dpois <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dpois",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape <- nimNumeric(J+1)
    priorRate <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape[j1] <<- model$getParam(marginalizedNodes[j1], 'shape') 
        priorRate[j1] <<- model$getParam(marginalizedNodes[j1], 'rate') 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + priorShape[j1] * log(priorRate[j1]) - (priorShape[j1] + y) * log(priorRate[j1] + 1) +
          lgamma(priorShape[j1] + y) - lgamma(priorShape[j1]) - lfactorial(y)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rgamma(1, shape = priorShape[j1] + y, rate = priorRate[j1] + 1))
      }
    }
  )
)


CRP_conjugate_dgamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape <- nimNumeric(J+1)
    priorRate <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape[j1] <<- model$getParam(marginalizedNodes[j1], 'shape') 
        priorRate[j1] <<- model$getParam(marginalizedNodes[j1], 'rate') 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + -0.5*log(2*pi) + priorShape[j1] * log(priorRate[j1]) - lgamma(priorShape[j1]) +
          lgamma(priorShape[j1] + 0.5) - (priorShape[j1] + 0.5)*log(priorRate[j1] + 0.5*(y-dataMean)^2)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataMean <- model$getParam(dataNodes[(i-1)*J+j1], 'mean')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rgamma(1, shape = priorShape[j1] + 0.5, rate = priorRate[j1] + (y-dataMean)^2/2))
      }
    }
  )
)


CRP_conjugate_dbeta_dbern <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dbern",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape1 <- nimNumeric(J+1)
    priorShape2 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape1[j1] <<- model$getParam(marginalizedNodes[j1], 'shape1') 
        priorShape2[j1] <<- model$getParam(marginalizedNodes[j1], 'shape2')
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + lgamma(priorShape1[j1]+y) + lgamma(priorShape2[j1]+1-y) - lgamma(priorShape1[j1]) - lgamma(priorShape2[j1]) - log(priorShape1[j1]+priorShape2[j1])
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rbeta(1, shape1=priorShape1[j1]+y, shape2=priorShape2[j1]+1-y)) 
      }
    }
  )
)

CRP_conjugate_dbeta_dbin <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dbin",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape1 <- nimNumeric(J+1)
    priorShape2 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape1[j1] <<- model$getParam(marginalizedNodes[j1], 'shape1') 
        priorShape2[j1] <<- model$getParam(marginalizedNodes[j1], 'shape2') 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        dataSize <- model$getParam(dataNodes[(i-1)*J+j1], 'size')
        out <- out + lgamma(priorShape1[j1]+priorShape2[j1]) + lgamma(priorShape1[j1]+y) + lgamma(priorShape2[j1]+dataSize-y) -
          lgamma(priorShape1[j1]) - lgamma(priorShape2[j1]) - lgamma(priorShape1[j1]+priorShape2[j1]+dataSize) +
          lfactorial(dataSize) - lfactorial(y) - lfactorial(dataSize-y)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataSize <- model$getParam(dataNodes[(i-1)*J+j1], 'size')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rbeta(1, shape1=priorShape1[j1]+y, shape2=priorShape2[j1]+dataSize-y))
      }
    }
  )
)


CRP_conjugate_dbeta_dnegbin <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dnegbin",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape1 <- nimNumeric(J+1)
    priorShape2 <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape1[j1] <<- model$getParam(marginalizedNodes[j1], 'shape1') 
        priorShape2[j1] <<- model$getParam(marginalizedNodes[j1], 'shape2') 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        dataSize <- model$getParam(dataNodes[(i-1)*J+j1], 'size')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + lgamma(priorShape1[j1]+priorShape2[j1]) + lgamma(priorShape1[j1]+dataSize) + lgamma(priorShape2[j1]+y) -
          lgamma(priorShape1[j1]) - lgamma(priorShape2[j1]) - lgamma(priorShape1[j1]+priorShape2[j1]+dataSize+y) +
          lfactorial(y+dataSize-1) - lfactorial(y) - lfactorial(dataSize-1)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataSize <- model$getParam(dataNodes[(i-1)*J+j1], 'size')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rbeta(1, shape1=priorShape1[j1]+dataSize, shape2=priorShape2[j1]+y))
      }
    }
  )
)

CRP_conjugate_dgamma_dexp <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dexp",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape <- nimNumeric(J+1)
    priorRate <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape[j1] <<- model$getParam(marginalizedNodes[j1], 'shape') 
        priorRate[j1] <<- model$getParam(marginalizedNodes[j1], 'rate')  
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + log(priorShape[j1]) + priorShape[j1]*log(priorRate[j1]) - (priorShape[j1]+1)*log(priorRate[j1]+y)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rgamma(1, shape=priorShape[j1]+1, rate=priorRate[j1]+y))
      }
    }
  )
)


CRP_conjugate_dgamma_dgamma <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dgamma",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape <- nimNumeric(J+1)
    priorRate <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape[j1] <<- model$getParam(marginalizedNodes[j1], 'shape') 
        priorRate[j1] <<- model$getParam(marginalizedNodes[j1], 'rate')   
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        datashape <- model$getParam(dataNodes[(i-1)*J+j1], 'shape')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + (datashape-1)*log(y) + priorShape[j1]*log(priorRate[j1]) + lgamma(datashape+priorShape[j1]) -
          lgamma(datashape) - lgamma(priorShape[j1]) -(datashape+priorShape[j1])*log(priorRate[j1]+y)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        datashape <- model$getParam(dataNodes[(i-1)*J+j1], 'shape')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rgamma(1, shape=datashape+priorShape[j1], rate=priorRate[j1]+y))
      }
    }
  )
)

CRP_conjugate_dgamma_dweib <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dweib",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape <- nimNumeric(J+1)
    priorRate <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape[j1] <<- model$getParam(marginalizedNodes[j1], 'shape') 
        priorRate[j1] <<- model$getParam(marginalizedNodes[j1], 'rate')   
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        dataShape <- model$getParam(dataNodes[(i-1)*J+j1], 'shape')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + log(dataShape) + (dataShape-1)*log(y) + priorShape[j1]*log(priorRate[j1]) +
          lgamma(priorShape[j1]+1) - lgamma(priorShape[j1]) - (priorShape[j1] + 1)*log(priorRate[j1] + y^dataShape)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataShape <- model$getParam(dataNodes[(i-1)*J+j1], 'shape')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rgamma(1, shape=1+priorShape[j1], rate=priorRate[j1]+y^dataShape))
      }
    }
  )
)

CRP_conjugate_dgamma_dinvgamma <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dinvgamma",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    priorShape <- nimNumeric(J+1)
    priorRate <- nimNumeric(J+1)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorShape[j1] <<- model$getParam(marginalizedNodes[j1], 'shape') 
        priorRate[j1] <<- model$getParam(marginalizedNodes[j1], 'rate')   
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        dataShape <- model$getParam(dataNodes[(i-1)*J+j1], 'shape')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        out <- out + -(dataShape+1)*log(y) + priorShape[j1]*log(priorRate[j1]) + lgamma(priorShape[j1] + dataShape) -
          lgamma(dataShape) - lgamma(priorShape[j1]) - (dataShape + priorShape[j1])*log(priorRate[j1] + 1/y)
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        dataShape <- model$getParam(dataNodes[(i-1)*J+j1], 'shape')
        y <- values(model, dataNodes[(i-1)*J+j1])[1]
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- c(rgamma(1, shape=dataShape+priorShape[j1], rate=priorRate[j1]+1/y))
      }
    }
  )
)

CRP_conjugate_ddirch_dmulti <- nimbleFunction(
  name = "CRP_conjugate_ddirch_dmulti",
  contains = CRP_helper,
  setup = function(model, marginalizedNodes, dataNodes, J, M) {
    d <- length(model[[marginalizedNodes[1]]])
    priorAlpha <- matrix(0, ncol=d, nrow=J)
  },
  methods = list(
    storeParams = function() {
      for(j1 in 1:J) {
        priorAlpha[j1, ] <<- model$getParam(marginalizedNodes[j1], 'alpha') 
      }
    },
    calculate_offset_coeff = function(i = integer(), j = integer()) {},
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      out <- 0
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        out <- out + lfactorial(d) - sum(lfactorial(y) + lgamma(priorAlpha[j1, ])) +
          lgamma(sum(priorAlpha[j1, ])) + sum(lgamma(priorAlpha[j1, ]+y)) - lgamma(sum(priorAlpha[j1, ]+y))
      }
      return(out)
    },
    sample = function(i = integer(), j = integer()) {
      for(j1 in 1:J) {
        y <- values(model, dataNodes[(i-1)*J+j1])
        values(model, marginalizedNodes[(j-1)*J+j1]) <<- rdirch(1, alpha = priorAlpha[j1, ]+y) 
      }
    }
  )
)


#' @rdname samplers
#' @export
sampler_CRP <- nimbleFunction(
  name = 'sampler_CRP',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    if(!is.null(control$printTruncation))
      printMessage <- control$printTruncation else printMessage <- TRUE

    targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    targetVar <- model$getVarNames(nodes = target)
    n <- length(targetElements)
    
    ## Find nodes indexed by the CRP node.
    if(is.null(control$clusterVarInfo)) {
        clusterVarInfo <- findClusterNodes(model, target)
    } else clusterVarInfo <- control$clusterVarInfo
    tildeVars <- clusterVarInfo$clusterVars

    ##  Various checks that model structure is consistent with our CRP sampler. 
    
    if(!is.null(clusterVarInfo$targetNonIndex))
      stop("sampler_CRP: Detected that the CRP variable is used in some way not as an index: ", clusterVarInfo$targetNonIndex, ". NIMBLE's CRP sampling not set up to handle this case.")
    
    if(any(clusterVarInfo$nTilde == 0)) {
      var <- which(clusterVarInfo$nTilde == 0)
      stop("sampler_CRP: Detected unusual indexing in ", deparse(clusterVarInfo$indexExpr[[var[1]]]),
           " . NIMBLE's CRP MCMC sampling is not designed for this case.")
    }
    
    if(is.null(tildeVars))
      stop('sampler_CRP: The model should have at least one cluster variable.\n')
    
    ## Cases like 'muTilde[xi[n-i+1]]'. sampler_CRP may be ok with this, but when we wrap the cluster node sampling
    ## to avoid sampling empty clusters, this kind of indexing will cause incorrect behavior.
    ## This case is trapped in findClusterNodes.
    
    allTildeNodes <- unlist(clusterVarInfo$clusterNodes)
    dataNodes <- model$getDependencies(target, stochOnly = TRUE, self = FALSE) # 'data' from the perspective of the clustering model
    stochDepsTildeNodes <- model$getDependencies(allTildeNodes, self = FALSE, stochOnly = TRUE)
    
    ## Make sure tildeNodes as determined from clustering actually are in model.
    if(!all(allTildeNodes %in% model$getNodeNames())) {
      missingNodes <- allTildeNodes[which(!allTildeNodes %in% model$getNodeNames())]
      stop("sampler_CRP: These cluster parameters are not nodes in the model: ", paste(missingNodes, collapse = ','))
    }
    
    ## Check that no other non-data nodes depend on cluster variables. 
    if(!identical(sort(dataNodes), sort(stochDepsTildeNodes)))
      stop("sampler_CRP: Only the variables being clustered can depend on the cluster parameters.")  

    ## Check that nodes in different clusters are distinct.
    ## E.g., this would not be the case if a user specified a joint prior on all cluster nodes
    ## Also check that clustering is done on stochastic nodes.
    for(varIdx in seq_along(clusterVarInfo$clusterVars)) {
        if(length(unique(clusterVarInfo$clusterNodes[[varIdx]])) != length(clusterVarInfo$clusterNodes[[varIdx]]))
            stop("sampler_CRP: cluster parameters in different clusters must be part of conditionally independent nodes.")
        if(any(model$isDeterm(clusterVarInfo$clusterNodes[[varIdx]])))
            stop("findClusterNodes: detected that deterministic nodes are being clustered. Please use the dCRP-distributed node to cluster stochastic nodes.")
    }
    
    ## Check that membership variable is independent of cluster nodes.
    ## Should be redundant with check that no other non-data nodes depend on cluster variables.
    if(target %in% stochDepsTildeNodes)
      stop("sampler_CRP: Cluster membership variable has to be independent of cluster parameters.")
    
    ## Check that cluster nodes are independent of membership variable
    ## (dataNodes are the dependents of target and should not contain cluster parameters).
    ## Should be redundant with check that no other non-data nodes depend on cluster variables.
    if(length(intersect(dataNodes, allTildeNodes)))
      stop("sampler_CRP: Cluster parameters have to be independent of cluster membership variable.")
    
    ## Check that observations are independent of each other.
    ## In non-conjugate case, this could potentially be relaxed within each cluster, provided we figure
    ## out correct ordering of dataNodes plus intermNodes in calculate().
    dataNodesIDs <- model$getDependencies(target, stochOnly = TRUE, self = FALSE, returnType = 'ids')
    sapply(dataNodes, function(x) {
      if(any(dataNodesIDs %in% model$getDependencies(x, self = FALSE, stochOnly = TRUE, returnType = 'ids')))
        stop("sampler_CRP: Variables being clustered must be conditionally independent. To model dependent variables being clustered jointly, you may use a multivariate distribution.")
    })

    ## Check that if same clusterNodes used twice in a declaration that they are used identically
    ## e.g., dnorm(thetaTilde[xi[i]], exp(thetaTilde[xi[i]])) is ok
    ## but dnorm(thetaTilde[xi[i]], exp(thetaTilde[xi[i]+1])) is not as can't properly add a new cluster
    ## (or account for when not to sample an empty cluster).
    if(length(unique(tildeVars)) != length(tildeVars)) 
        for(idx1 in 2:length(tildeVars))
            for(idx2 in 1:(length(tildeVars)-1))
                if(!identical(clusterVarInfo$clusterNodes[[idx1]], clusterVarInfo$clusterNodes[[idx2]]) &&
                   any(clusterVarInfo$clusterNodes[[idx1]] %in% clusterVarInfo$clusterNodes[[idx2]]))
                    stop("sampler_CRP: Inconsistent indexing in or inconsistent dependencies of ",
                         deparse(clusterVarInfo$indexExpr[[idx1]]), " and ",
                         deparse(clusterVarInfo$indexExpr[[idx2]]), ".")
    
    nData <- length(dataNodes)

    ## Check that no use of multiple clustering variables, such as 'thetaTilde[xi[i], eta[j]]'.
    ## It's likely that if we set the non-conjugate sampler and turn off wrapping omits sampling of empty clusters
    ## (which is not set up correctly for this case), that the existing code would give correct sampling.
    if(any(clusterVarInfo$multipleStochIndexes))
        stop("sampler_CRP: Detected use of multiple stochastic indexes of a variable: ", deparse(clusterVarInfo$indexExpr[[1]]), ". NIMBLE's CRP sampling is not yet set up to handle this case. Please contact the NIMBLE development team if you are interested in this functionality.")

    ## Check there is at least one or more "observation" per random index.
    ## Note that cases like mu[xi[i],xi[j]] are being trapped in findClusterNodes().
    if(n > nData)
       stop("sampler_CRP: At least one variable has to be clustered for each cluster membership ID.")
    
   
    p <- length(tildeVars) 
    nTilde <- clusterVarInfo$nTilde / clusterVarInfo$numNodesPerCluster  ## implied number of potential clusters
    if(length(unique(nTilde)) != 1)
        stop('sampler_CRP: In a model with multiple cluster parameters, the number of those parameters must all be the same.\n')
    min_nTilde <- nTilde[1]
    if(min_nTilde < n)
      warning('sampler_CRP: The number of clusters based on the cluster parameters is less than the number of potential clusters. The MCMC is not strictly valid if it ever proposes more components than cluster parameters exist; NIMBLE will warn you if this occurs.\n')
    
    ## Determine if concentration parameter is fixed or random (code similar to the one in sampleDPmeasure function).
    ## This is used in truncated case to tell user if model is proper or not.
    fixedConc <- TRUE
    parentNodesTarget <- getParentNodes(target, model, stochOnly = TRUE)
    if(length(parentNodesTarget)) {
      fixedConc <- FALSE
    }

    ## Here we try to set up some data structures that allow us to do observation-specific
    ## computation, to save us from computing for all observations when a single cluster membership is being proposed.
    ## At the moment, this is the only way we can easily handle dependencies for multiple node elements in a
    ## 'vectorized' way.
    nObsPerClusID <- nData / n # equal to one in standard CRP model (former J)
    dataNodes <- rep(targetElements[1], nData) ## this serves as dummy nodes that may be replaced below
    ## needs to be legitimate nodes because run code sets up calculate even if if() would never cause it to be used
    for(i in seq_len(n)) { # dataNodes are always needed so only create them before creating  intermNodes
      stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self = FALSE)
      dataNodes[((i-1)*nObsPerClusID + 1) : (i*nObsPerClusID)] <- stochDeps
    }
    nIntermClusNodesPerClusID <- length(model$getDependencies(targetElements[1], determOnly = TRUE))  #nInterm
    intermNodes <- dataNodes # initialize nodes in case not used (needed for compilation to go through)
    if(nIntermClusNodesPerClusID > 0) {
      intermNodes <- rep(as.character(NA), nIntermClusNodesPerClusID * n)
      for(i in seq_len(n)) {
        detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
        intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)] <- detDeps
      }
    }
    
    helperFunctions <- nimbleFunctionList(CRP_helper)
        
    ## use conjugacy to determine which helper functions to use
    if(control$checkConjugacy) {
        conjugacyResult <- checkCRPconjugacy(model, target) 
    } else conjugacyResult <- NULL
    
    if(is.null(conjugacyResult)) {
      sampler <- 'CRP_nonconjugate'
    } else 
      sampler <- switch(conjugacyResult,
                        conjugate_dnorm_dnorm = 'CRP_conjugate_dnorm_dnorm',
                        conjugate_dmnorm_dmnorm = 'CRP_conjugate_dmnorm_dmnorm',
                        conjugate_dnorm_dnorm_nonidentity = 'CRP_conjugate_dnorm_dnorm_nonidentity',
                        conjugate_dinvgamma_dnorm = 'CRP_conjugate_dinvgamma_dnorm',
                        conjugate_dinvwish_dmnorm = 'CRP_conjugate_dinvwish_dmnorm',
                        conjugate_dwish_dmnorm = 'CRP_conjugate_dwish_dmnorm',
                        conjugate_dnorm_invgamma_dnorm = 'CRP_conjugate_dnorm_invgamma_dnorm',
                        conjugate_dnorm_gamma_dnorm = 'CRP_conjugate_dnorm_gamma_dnorm',
                        conjugate_dmnorm_invwish_dmnorm = 'CRP_conjugate_dmnorm_invwish_dmnorm',
                        conjugate_dmnorm_wish_dmnorm = 'CRP_conjugate_dmnorm_wish_dmnorm',
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


    clusterIDs <- unique(clusterVarInfo$clusterIDs[[1]])
    nClusters <- length(clusterVarInfo$clusterIDs)
    if(nClusters > 1) {
        ## Check that each set of tildeNodes indicate same number of clusters.        
        sapply(2:nClusters, function(i) {
            if(!identical(clusterIDs, unique(clusterVarInfo$clusterIDs[[i]])))
                stop("sampler_CRP: differing number of clusters indicated by ", paste0(clusterVarInfo$clusterNodes[[1]], collapse = ', '), " and ", paste0(clusterVarInfo$clusterNodes[[i]], collapse = ', '), ".")
        })
    }

    nClusNodesPerClusID <- sum(clusterVarInfo$numNodesPerCluster)

    ## Determine correct order of clusterNodes, including any intermediate nodes
    ## standing between different clusterNodes. This set of nodes are the
    ## 'marginalized' nodes.
    ## Also check independence of cluster parameters across clusters.
    ## This block seems to increase buildMCMC time by about 30%.
    allNodes <- unlist(clusterVarInfo$clusterNodes)
    marginalizedNodes <- allNodes
    dataIntermNodes <- c(dataNodes, intermNodes)
    ids <- unlist(clusterVarInfo$clusterIDs)
    if(nClusNodesPerClusID > 1) {
        for(id in seq_along(clusterIDs)) {
            origNodes <- allNodes[ids == id]
            nodes <- model$getDependencies(origNodes)
            nodes <- nodes[!nodes %in% dataIntermNodes]        
            if(id == min(ids)) {
                nClusNodesPerClusID <- length(nodes)  # now includes determ intermediates between cluster nodes
                marginalizedNodes <- rep('', nClusNodesPerClusID * min_nTilde)
            }
            if(length(nodes) != nClusNodesPerClusID)
                stop("sampler_CRP: detected differing number of cluster parameters across clusters or dependence of parameters across clusters.") # unequal number could occur because of dependence across clusters
            marginalizedNodes[((id-1)*nClusNodesPerClusID+1):(id*nClusNodesPerClusID)] <- nodes
            if(any(nodes %in% allNodes[!allNodes %in% origNodes]))
                stop("sampler_CRP: cluster parameters must be independent across clusters.")
        }
        
    } else {
        marginalizedNodes <- clusterVarInfo$clusterNodes[[1]]
        deps <- unique(unlist(lapply(marginalizedNodes, function(x)
            model$getDependencies(x, self = FALSE, stochOnly = TRUE))))
        if(any(marginalizedNodes %in% deps))
            stop("sampler_CRP: cluster parameters must be independent across clusters.")
        if(!identical(sort(deps), sort(dataNodes)))
            warning("sampler_CRP: dependencies of cluster parameters include unexpected nodes: ",
                    paste0(deps[!deps %in% dataNodes], collapse = ', '))
    }
    
    identityLink <- TRUE

    if(p == 2 && sampler %in% c("CRP_conjugate_dnorm_invgamma_dnorm", "CRP_conjugate_dnorm_gamma_dnorm", 
                                "CRP_conjugate_dmnorm_invwish_dmnorm", "CRP_conjugate_dmnorm_wish_dmnorm")) {
        if(sampler == "CRP_conjugate_dnorm_invgamma_dnorm") {
            dist1 <- 'dnorm'
            dist2 <- 'dinvgamma'
        }
        if(sampler == "CRP_conjugate_dnorm_gamma_dnorm") {
            dist1 <- 'dnorm'
            dist2 <- 'dgamma'
        }
        if(sampler == "CRP_conjugate_dmnorm_invwish_dmnorm") {
            dist1 <- 'dmnorm'
            dist2 <- 'dinvwish'
        }
        if(sampler == "CRP_conjugate_dmnorm_wish_dmnorm") {
            dist1 <- 'dmnorm'
            dist2 <- 'dwish'
        }
        for(i in seq_along(tildeVars)) {
            if(model$getDistribution(clusterVarInfo$clusterNodes[[i]][1]) == dist1) {
                marginalizedNodes1 <- clusterVarInfo$clusterNodes[[i]]
            } 
            if(model$getDistribution(clusterVarInfo$clusterNodes[[i]][1]) == dist2) {
                marginalizedNodes2 <- clusterVarInfo$clusterNodes[[i]]
            }
        }  
      helperFunctions[[1]] <- eval(as.name(sampler))(model, marginalizedNodes1, marginalizedNodes2, dataNodes, nObsPerClusID, nClusNodesPerClusID)
      calcNodes <- model$getDependencies(c(target, marginalizedNodes1, marginalizedNodes2))
    } else {
      calcNodes <- model$getDependencies(c(target, marginalizedNodes))
      if(sampler == "CRP_conjugate_dnorm_dnorm_nonidentity") {
          identityLink <- FALSE
          helperFunctions[[1]] <- eval(as.name(sampler))(model, marginalizedNodes, dataNodes, intermNodes, nIntermClusNodesPerClusID, nObsPerClusID, nClusNodesPerClusID)
      } else {
          if(sampler != "CRP_nonconjugate" && nObsPerClusID != nClusNodesPerClusID)
              ## Our conjugate samplers use J (nObsPerClusID) instead of M (nClusNodesPerClusID) because they assume J=M.
              ## We shouldn't get to this point (based on conjugacy checking), but catch this if we do.
              stop("Number of observations per group does not equal number of cluster parameters per group. NIMBLE's CRP sampling is not set up to handle this except for the non-conjugate sampler.")
          helperFunctions[[1]] <- eval(as.name(sampler))(model, marginalizedNodes, dataNodes, nObsPerClusID, nClusNodesPerClusID)
      }
    }

    curLogProb <- numeric(n)
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
      if(printMessage) {
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
            model[[target]][i] <<- xiUniques[j] # <<-
            if(nIntermClusNodesPerClusID > 0) {
              model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
            }
            curLogProb[iprob] <<- log(xiCounts[xiUniques[j]]) +
                model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])
            reorderXiUniques[iprob] <- xiUniques[j]
            iprob <- iprob + 1
          }
        }
        
        ## Second, compute probability of sampling a new cluster, here, new cluster is the current cluster!
        model[[target]][i] <<- xi[i] # <<- label of new component
        if(sampler == 'CRP_nonconjugate'){ # simulate tildeVars[xi[i]] # do this everytime there is a singleton so we ensure this comes always from the prior
          helperFunctions[[1]]$sample(i, model[[target]][i])
          if(nIntermClusNodesPerClusID > 0) {
            model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
          }
          model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])       
        }
        if(!identityLink) 
            helperFunctions[[1]]$calculate_offset_coeff(i, model[[target]][i])
        curLogProb[k] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i) # <<- probability of sampling a new label, only k components because xi_i is a singleton
        
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
          if(nIntermClusNodesPerClusID > 0) {
            model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
          }
          curLogProb[j] <<- log(xiCounts[xiUniques[j]]) +
              model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])       
        }
        ## Second, compute probability of sampling a new cluster depending on the value of kNew.       
        if(kNew == 0) { # no new cluster can be created 
          curLogProb[k+1] <<- log(0)  # <<- k+1 <= n always because k==n requires all singletons, handled above
        } else { # a new cluster can be created
          model[[target]][i] <<- kNew 
          if(sampler == 'CRP_nonconjugate'){
            helperFunctions[[1]]$sample(i, model[[target]][i])
            if(nIntermClusNodesPerClusID > 0) {
              model$calculate(intermNodes[((i-1)*nIntermClusNodesPerClusID+1):(i*nIntermClusNodesPerClusID)]) 
            }
            model$calculate(dataNodes[((i-1)*nObsPerClusID+1):(i*nObsPerClusID)])       
          }
          if(!identityLink) 
            helperFunctions[[1]]$calculate_offset_coeff(i, model[[target]][i])
          curLogProb[k+1] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i) # <<- probability of sampling a new label
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
        ## Reset to previous marginalized node value; we choose to store information on what elements to be restored in sample()
        ## but an alternative would be to have i=0 determine reset and pass j=kNew here.
        if(sampler == 'CRP_nonconjugate')   
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
  )
)

findClusterNodes <- function(model, target) {
  ## Determine which model nodes are the cluster parameters by processing expressions to look
  ## for what is indexed by the dCRP clusterID nodes. This also determine which clusterID
  ## each cluster parameter is associated with.
  targetVar <- model$getVarNames(nodes = target)
  targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
  deps <- model$getDependencies(target, self = FALSE, returnType = 'ids')
  declIDs <- model$modelDef$maps$graphID_2_declID[deps] ## declaration IDs of the nodeIDs
  uniqueIDs <- unique(declIDs)
  depsByDecl <- lapply(uniqueIDs, function(x) deps[which(x == declIDs)])

  ## Find one example dependency per BUGS declaration for more efficient processing
  exampleDeps <- model$modelDef$maps$graphID_2_nodeName[sapply(depsByDecl, `[`, 1)]
    
  ## Once we find the cluster parameter variables below, we want to evaluate the cluster membership
  ## values (e.g., xi[1],...,xi[n]) for all possible values they could take, this will
  ## allow us to determine all possible cluster nodes in the model (though some may
  ## not actually be specified in the model, if there is truncation).
  ## Therefore, set up an evaluation environment in which (xi[1],...,xi[n]) = (1,2,...,n)
  ## first try was: e[[targetVar]] <- seq_along(targetElements)
  ## However in first try, that wouldn't handle xi[3:10] ~ dCRP(), but next construction does.
  e <- list()
  idxExpr <- model$getDeclInfo(target)[[1]]$indexExpr[[1]]
  eval(substitute(`<-`(`[`(e$VAR, IDX), seq_along(targetElements)), list(VAR = targetVar, IDX = idxExpr)))
  ## For cases of cross clustering (e.g., mu[xi[i],eta[j]]) we need the other dcrp node(s)
  nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  dists <- model$getDistribution(nodes)
  if(length(dists == 'dCRP') > 1) { 
      dcrpNodes <- nodes[dists == 'dCRP' & nodes != target]
      for(i in seq_along(dcrpNodes)) {
          dcrpElements <- model$expandNodeNames(dcrpNodes[i], returnScalarComponents = TRUE)
          dcrpVar <- model$getVarNames(nodes = dcrpNodes[i])
          idxExpr <- model$getDeclInfo(dcrpNodes[i])[[1]]$indexExpr[[1]]
          eval(substitute(`<-`(`[`(e$VAR, IDX), seq_along(dcrpElements)), list(VAR = dcrpVar, IDX = idxExpr)))

      }
  }
  clusterNodes <- indexExpr <- clusterIDs <- list()
  clusterVars <- indexPosition <- numIndexes <- targetIsIndex <- targetIndexedByFunction <-
      loopIndex <- NULL
  varIdx <- 0

  targetNonIndex <- NULL
  multipleStochIndexes <- NULL

  modelVars <- model$getVarNames()
  modelVars <- modelVars[!modelVars == targetVar]

  ## Process model declaration expressions to find stochastic indexing and the indexed variable.
  for(idx in seq_along(exampleDeps)) {
    ## Pull out expressions, either as RHS of deterministic or parameters of stochastic
    fullExpr <- cc_getNodesInExpr(model$getValueExpr(exampleDeps[idx]))
    for(j in seq_along(fullExpr)) {
      subExpr <- parse(text = fullExpr[j])[[1]]  # individual parameter of stochastic or RHS of deterministic
      len <- length(subExpr)
      ## Look for target variable within expression, but only when used within index
      if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[' &&
        sum(all.vars(subExpr) == targetVar) && subExpr[[2]] != targetVar) {
        varIdx <- varIdx + 1
        multipleStochIndexes <- c(multipleStochIndexes, FALSE)
          
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
            ## We will need to relax this when allow crossed clustering.
            if(sum(all.vars(subExpr[[k]]) %in% modelVars)) { ## cases like mu[xi[i],eta[j]]
                ## We are adding support for this case.
                ## warning("findClusterNodes: multiple indexing variables in '", deparse(subExpr),
                ##          "'. NIMBLE's CRP MCMC sampling not designed for this situation.")
                multipleStochIndexes[varIdx] <- TRUE
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

        if(targetIndexedByFunction[varIdx] && ncol(unrolledIndices) > 1)  ## Now that we allow cluster parameters with multiple indexes, this is very hard to handle in terms of identifying what column of unrolledIndices to use for sorting clusterNodes.
            stop("findClusterNodes: Detected that a cluster parameter is indexed by a function such as 'mu[xi[n-i+1]]' rather than simple indexing such as 'mu[xi[i]]'. NIMBLE's CRP MCMC sampling not designed for this case.")
        loopIndexes <- unlist(sapply(declInfo$symbolicParentNodes,
                                    function(x) {
                                        if(length(x) >= 3 && x[[1]] == '[' &&
                                           x[[2]] == targetVar) return(deparse(x[[3]]))
                                        else return(NULL) }))
        if(length(loopIndexes) != 1)
            stop("findClusterNodes: found cluster membership parameters that use different indexing variables; NIMBLE's CRP sampling not designed for this case.")
        ## Note not clear when NULL would be the result...
        loopIndex[varIdx] <- loopIndexes

        ## Determine potential cluster nodes by substituting all possible clusterID values into the indexing expression. 
        n <- nrow(unrolledIndices)
        if(n > 0 && loopIndex[varIdx] %in% dimnames(unrolledIndices)[[2]]) {  # catch cases like use of xi[2] rather than xi[i]
            ## Order so that loop over index of cluster ID in order of cluster ID so that
            ## clusterNodes will be grouped in chunks of unique cluster IDs for correct
            ## sampling of new clusters when have multiple obs per cluster.
            ord <- order(unrolledIndices[ , loopIndex[varIdx]])
            unrolledIndices <- unrolledIndices[ord, , drop = FALSE]
            clusterIDs[[varIdx]] <- unrolledIndices[ , loopIndex[varIdx]]

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
        } else {
            clusterNodes[[varIdx]] <- character(0)
            clusterIDs[[varIdx]] <- numeric(0)
        }
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
    if(nTilde[varIdx]) {
        if(any(is.na(clusterNodes[[varIdx]])))  
            stop("findClusterNodes: fewer cluster IDs in ", target, " than elements being clustered.")
        
        ## Handle cases where indexing of variables in dynamic indexing does not correspond to actual
        ## stochastic model nodes.
        if(any(!clusterNodes[[varIdx]] %in% modelNodes)) {
            clusterNodes[[varIdx]] <- lapply(clusterNodes[[varIdx]], function(x) model$expandNodeNames(x))
            clusterIDs[[varIdx]] <- rep(clusterIDs[[varIdx]], times = sapply(clusterNodes[[varIdx]], length))
            clusterNodes[[varIdx]] <- unlist(clusterNodes[[varIdx]])
        }
        ## Now remove duplicates when indexed variables correspond to same model node,
        ## but only for duplicates within a cluster.
        groups <- split(clusterNodes[[varIdx]], clusterIDs[[varIdx]])
        dups <- unlist(lapply(groups, duplicated))
        clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][!dups]
        clusterIDs[[varIdx]] <- clusterIDs[[varIdx]][!dups]

        ## Formerly we were checking that we had a contiguous set of cluster nodes
        ## starting with the first one, but for clusterNodes with more than one index and
        ## truncation this is hard to do, so just fall back to returning the clusterNodes
        ## that are actually part of the model.
        validNodes <- clusterNodes[[varIdx]] %in% modelNodes
        
        if(!all(validNodes)) {  # i.e., truncated representation
            clusterNodes[[varIdx]] <- clusterNodes[[varIdx]][validNodes]
            clusterIDs[[varIdx]] <- clusterIDs[[varIdx]][validNodes]
        }
    }
  }

  nTilde <- sapply(clusterNodes, length)
  numNodesPerCluster <- sapply(clusterIDs, function(x) {
      tbl <- table(x)
      num <- unique(tbl)
      if(length(num) > 1) stop("findClusterNodes: detected differing numbers of nodes (i.e., parameters) per cluster. NIMBLE's CRP sampling not designed for this case.")
      return(num)})
      
  return(list(clusterNodes = clusterNodes, clusterVars = clusterVars, nTilde = nTilde,
              numNodesPerCluster = numNodesPerCluster, clusterIDs = clusterIDs, loopIndex = loopIndex,
              targetIsIndex = targetIsIndex, indexPosition = indexPosition, indexExpr = indexExpr,
              numIndexes = numIndexes, targetIndexedByFunction = targetIndexedByFunction,
              targetNonIndex = targetNonIndex, multipleStochIndexes = multipleStochIndexes))
}


checkCRPconjugacy <- function(model, target) {
    ## Checks if can use conjugacy in drawing new components for dCRP node updating.
    ## Should detect various univariate and multivariate cases.
    ## We currently handle only a limited dnorm-dnorm non-identity relationship,
    ## e.g., b0[xi[i]] + b1*x[i] or b0 + b1[xi[i]]*x[i], but not b0[xi[i]] + b1[xi[i]]*x[i]
    ## We plan to add dpois-dgamma non-identity relationship.
    
    conjugate <- FALSE 

    targetAsScalars <- model$expandNodeNames(target, returnScalarComponents=TRUE) 
    targetElementExample <- targetAsScalars[1]
    n <- length(targetAsScalars)

    clusterVarInfo <- findClusterNodes(model, target)
 
    ## Check conjugacy for one cluster node (for efficiency reasons) and then make sure all cluster nodes are IID.
    ## Since we only allow one clusterVar, shouldn't need to worry that depNodes for difference clusters are
    ## from different declarations (e.g.  y[1] ~ dnorm(thetaTilde[xi[1]],1) and y[2] ~ dt(thetaTilde[xi[2]],1,1).
    if(length(clusterVarInfo$clusterVars) == 1) {  ## for now avoid case of mixing over multiple parameters, but allow dnorm_dinvgamma below
        clusterNodes <- clusterVarInfo$clusterNodes[[1]]  # e.g., 'thetatilde[1]',...,
        clusterIDs <- clusterVarInfo$clusterIDs[[1]]
        clusterNodesFirst <- clusterNodes[clusterIDs == 1] # need to check for all nodes in a cluster
        ## Currently we only handle offsets and coeffs for dnorm case;
        ## will add Pois-gamma and possibly MVN cases.
        identityLink <- TRUE
        conjugacy <- model$checkConjugacy(clusterNodesFirst, restrictLink = 'identity')
        if(!(length(conjugacy) == length(clusterNodesFirst)) && all(model$getDistribution(clusterNodesFirst) == 'dnorm')) {
            identityLink <- FALSE
            conjugacy <- model$checkConjugacy(clusterNodesFirst)  ## check non-identity link too
        }
        if(length(conjugacy) == length(clusterNodesFirst) && length(unique(sapply(conjugacy, '[[', 'type'))) == 1) {
            ## All conjugate and all same conjugacy type.
            conjugacyType <- paste0(conjugacy[[1]]$type, '_', sub('dep_', '', names(conjugacy[[1]]$control)))
            if(!identityLink)
                conjugacyType <- paste0(conjugacyType, '_nonidentity')
            conjugate <- TRUE
            ## Check that dependent nodes ('observations') from same declaration.
            ## This should ensure they have same distribution and parameters are being
            ## clustered in same way, but also allows other parameters to vary, e.g.,
            ## y[i] ~ dnorm(mu[xi[i]], s2[i])
            depNodes <- model$getDependencies(clusterNodesFirst, stochOnly = TRUE, self=FALSE)
            if(length(unique(model$getDeclID(depNodes))) != 1)  ## make sure all dependent nodes from same declaration (i.e., exchangeable)
                conjugate <- FALSE
            ## We need each data node to have corresponding cluster parameter.
            if(length(depNodes) / n != length(clusterNodesFirst))
                conjugate <- FALSE

            ## Check that cluster nodes are IID (across clusters), since for efficiency we only
            ## check the conjugacy for the first cluster above.
            ## Extended to work with models with more than one observation per cluster ID.
            splitNodes <- split(clusterNodes, clusterIDs)
            valueExprs <- lapply(splitNodes, function(x) {
                out <- sapply(x, model$getValueExpr)
                names(out) <- NULL
                out
            })
            if(length(unique(valueExprs)) != 1) 
                conjugate <- FALSE
        }
    }
    ## check for dnorm_dinvgamma conjugacy
    if(length(clusterVarInfo$clusterVars) == 2 &&
      checkNormalInvGammaConjugacy(model, clusterVarInfo, n, 'dgamma')) {
        conjugate <- TRUE
        conjugacyType <- "conjugate_dnorm_gamma_dnorm"
    } else if(length(clusterVarInfo$clusterVars) == 2 &&
      checkNormalInvGammaConjugacy(model, clusterVarInfo, n, 'dinvgamma')) {
        conjugate <- TRUE
        conjugacyType <- "conjugate_dnorm_invgamma_dnorm"
    } else if(length(clusterVarInfo$clusterVars) == 2 &&
      checkNormalInvWishartConjugacy(model, clusterVarInfo, n, 'dwish')) {
        conjugate <- TRUE
        conjugacyType <- "conjugate_dmnorm_wish_dmnorm"
    } else if(length(clusterVarInfo$clusterVars) == 2 &&
      checkNormalInvWishartConjugacy(model, clusterVarInfo, n, 'dinvwish')) {
        conjugate <- TRUE
        conjugacyType <- "conjugate_dmnorm_invwish_dmnorm"
    }
    if(conjugate) return(conjugacyType) else return(NULL)
}


checkNormalInvGammaConjugacy <- function(model, clusterVarInfo, n, gammaDist = 'dinvgamma') {
    ## This function can also check dnorm_gamma when 'gammaDist' is 'dgamma'.
    if(length(clusterVarInfo$clusterVars) != 2)
        stop("checkNormalInvGammaConjugacy: requires two cluster variables.")
    conjugate <- FALSE

    varParam <- 'var'
    if(gammaDist == 'dgamma')
        varParam <- 'tau'
    
    clusterNodes1 <- clusterVarInfo$clusterNodes[[1]]
    clusterNodes2 <- clusterVarInfo$clusterNodes[[2]]

    if(!any(model$isDeterm(c(clusterNodes1, clusterNodes2))) &&
       length(clusterNodes1) == length(clusterNodes2) &&
       identical(clusterVarInfo$clusterIDs[[1]], clusterVarInfo$clusterIDs[[2]])) {
        dists <- c(model$getDistribution(clusterNodes1[1]), model$getDistribution(clusterNodes2[1]))
        if(dists[1] == gammaDist && dists[2] ==  "dnorm") {  ## put in order so dnorm node is first
            dists <- c("dnorm", gammaDist)
            tmp <- clusterNodes1; clusterNodes1 <- clusterNodes2; clusterNodes2 <- tmp
        }

        clusterIDs <- clusterVarInfo$clusterIDs[[1]]
        
        ## Check conjugacy for example nodes.
        exampleNodes1 <- clusterNodes1[clusterIDs == 1]
        exampleNodes2 <- clusterNodes2[clusterIDs == 1]
        if(dists[1] == "dnorm" && dists[2] == gammaDist) {
            conjugacy_dnorm <- model$checkConjugacy(exampleNodes1, restrictLink = 'identity')
            conjugacy_dinvgamma <- model$checkConjugacy(exampleNodes2)
            if(length(conjugacy_dnorm) == length(exampleNodes1) &&
               length(conjugacy_dinvgamma) == length(exampleNodes2) &&
               sum(sapply(conjugacy_dnorm, '[[', 'type') == 'conjugate_dnorm') == length(exampleNodes1) &&
               sum(sapply(conjugacy_dinvgamma, '[[', 'type') == paste0('conjugate_', gammaDist)) == length(exampleNodes2) &&
               all(sapply(seq_along(conjugacy_dinvgamma), function(idx)
                   sum(conjugacy_dinvgamma[[idx]]$control$dep_dnorm == exampleNodes1[idx]) == 1)))
                conjugate <- TRUE
        }
        if(conjugate) {
            ## Check that cluster nodes are IID (across clusters), since for efficiency above
            ## we only check the fist cluster. Can have non-IID within a cluster.
            if(any(model$getDistribution(clusterNodes1) != "dnorm"))
                conjugate <- FALSE

            splitNodes1 <- split(clusterNodes1, clusterIDs)
            splitNodes2 <- split(clusterNodes2, clusterIDs)
            ## Check that means of cluster mean nodes are same across clusters.
            meanExprs <- lapply(splitNodes1, function(x) {
                out <- sapply(x, function(z) model$getParamExpr(z, 'mean'))
                names(out) <- NULL
                out
            })
            if(length(unique(meanExprs)) != 1)
                conjugate <- FALSE
            ## Check that variance for cluster mean nodes are same except for dependence on variance.
            varExprs <- lapply(seq_along(splitNodes1), function(idx) {
                exprs <- sapply(splitNodes1[[idx]], function(z) model$getParamExpr(z, varParam))
                exprs <- sapply(exprs, function(expr) cc_expandDetermNodesInExpr(model, expr))
                names(exprs) <- NULL
                for(i in seq_along(exprs)) {
                    varText <- deparse(exprs[[i]])
                    if(length(grep(splitNodes2[[idx]][i], varText, fixed = TRUE)))  ## remove clusterNodes2[i] from expression so var expressions will be the same
                        exprs[[i]] <- parse(text = gsub(splitNodes2[[idx]][i],
                                                        "a", varText, fixed = TRUE))[[1]]
                }
                exprs
            })
            if(length(unique(varExprs)) != 1) 
                conjugate <- FALSE
            
            ## Check that cluster variance nodes are IID
            valueExprs <- lapply(splitNodes2, function(x) {
                out <- sapply(x, model$getValueExpr)
                names(out) <- NULL
                out
            })
            if(length(unique(valueExprs)) != 1)
                conjugate <- FALSE
            
            ## Check that dependent nodes ('observations') from same declaration.
            ## This should ensure they have same distribution and parameters are being
            ## clustered in same way.
            depNodes <- model$getDependencies(exampleNodes1, stochOnly = TRUE, self = FALSE)
            if(length(unique(model$getDeclID(depNodes))) != 1)
                conjugate <- FALSE
            ## We are not set up to have multiple data nodes depend on a single normal-invgamma distribution
            if(length(exampleNodes1) != length(depNodes) / n)
                conjugate <- FALSE
        }
    }
    return(conjugate)
}
    
checkNormalInvWishartConjugacy <- function(model, clusterVarInfo, n, wishartDist = 'dinvwish') {
    ## This function can also check dmnorm_wish when 'wishartDist' is 'dwish'.
    if(length(clusterVarInfo$clusterVars) != 2)
        stop("checkNormalInvWishartConjugacy: requires two cluster variables.")
    conjugate <- FALSE
    
    covParam <- 'cov'
    if(wishartDist == 'dwish')
        covParam <- 'prec'

    clusterNodes1 <- clusterVarInfo$clusterNodes[[1]]
    clusterNodes2 <- clusterVarInfo$clusterNodes[[2]]

    if(!any(model$isDeterm(c(clusterNodes1, clusterNodes2))) &&
       length(clusterNodes1) == length(clusterNodes2) &&
       identical(clusterVarInfo$clusterIDs[[1]], clusterVarInfo$clusterIDs[[2]])) {
        dists <- c(model$getDistribution(clusterNodes1[1]), model$getDistribution(clusterNodes2[1]))
        if(dists[1] == wishartDist && dists[2] ==  "dmnorm") {  ## put in order so dmnorm node is first
            dists <- c("dmnorm", wishartDist)
            tmp <- clusterNodes1; clusterNodes1 <- clusterNodes2; clusterNodes2 <- tmp
        }

        clusterIDs <- clusterVarInfo$clusterIDs[[1]]
        
        ## Check conjugacy for example nodes.
        exampleNodes1 <- clusterNodes1[clusterIDs == 1]
        exampleNodes2 <- clusterNodes2[clusterIDs == 1]

        if(dists[1] == "dmnorm" && dists[2] == wishartDist) {
            conjugacy_dmnorm <- model$checkConjugacy(exampleNodes1, restrictLink = 'identity')
            conjugacy_dinvwish <- model$checkConjugacy(exampleNodes2)
            if(length(conjugacy_dmnorm) == length(exampleNodes1) &&
               length(conjugacy_dinvwish) == length(exampleNodes2) &&
               sum(sapply(conjugacy_dmnorm, '[[', 'type') == 'conjugate_dmnorm') == length(exampleNodes1) &&
               sum(sapply(conjugacy_dinvwish, '[[', 'type') == paste0('conjugate_', wishartDist)) == length(exampleNodes2) &&
               all(sapply(seq_along(conjugacy_dinvwish), function(idx)
                   sum(conjugacy_dinvwish[[idx]]$control$dep_dmnorm == exampleNodes1[idx]) == 1)))
                conjugate <- TRUE
        }
        if(conjugate) {  
            ## Check that cluster nodes are IID (across clusters), since for efficiency above
            ## we only check the fist cluster. Can have non-IID within a cluster.
            if(any(model$getDistribution(clusterNodes1) != "dmnorm"))
                conjugate <- FALSE
            
            splitNodes1 <- split(clusterNodes1, clusterIDs)
            splitNodes2 <- split(clusterNodes2, clusterIDs)
            ## Check that means of cluster mean nodes are same across clusters.
            meanExprs <- lapply(splitNodes1, function(x) {
                out <- sapply(x, function(z) model$getParamExpr(z, 'mean'))
                names(out) <- NULL
                out
            })
            if(length(unique(meanExprs)) != 1)
                conjugate <- FALSE
            ## Check that variance for cluster mean nodes are same except for dependence on variance.
            varExprs <- lapply(seq_along(splitNodes1), function(idx) {
                exprs <- sapply(splitNodes1[[idx]], function(z) model$getParamExpr(z, covParam))
                exprs <- sapply(exprs, function(expr) cc_expandDetermNodesInExpr(model, expr))
                names(exprs) <- NULL
                for(i in seq_along(exprs)) {
                    varText <- deparse(exprs[[i]])
                    if(length(grep(splitNodes2[[idx]][i], varText, fixed = TRUE)))  ## remove clusterNodes2[i] from expression so var expressions will be the same
                        exprs[[i]] <- parse(text = gsub(splitNodes2[[idx]][i],
                                                        "a", varText, fixed = TRUE))[[1]]
                }
                exprs
            })
            if(length(unique(varExprs)) != 1) 
                conjugate <- FALSE
            
            ## Check that cluster variance nodes are IID
            valueExprs <- lapply(splitNodes2, function(x) {
                out <- sapply(x, model$getValueExpr)
                names(out) <- NULL
                out
            })
            if(length(unique(valueExprs)) != 1)
                conjugate <- FALSE
            
            ## Check that dependent nodes ('observations') from same declaration.
            ## This should ensure they have same distribution and parameters are being
            ## clustered in same way.
            depNodes <- model$getDependencies(clusterNodes1, stochOnly = TRUE, self = FALSE)
            if(length(unique(model$getDeclID(depNodes))) != 1)
                conjugate <- FALSE
            ## We are not set up to have multiple data nodes depend on a single normal-invwish distribution
            if(length(exampleNodes1) != length(depNodes) / n)
                conjugate <- FALSE
        }
    }
    return(conjugate)
}

sampler_CRP_cluster_wrapper <- nimbleFunction(
    name = "CRP_cluster_wrapper", 
    contains = sampler_BASE,
    setup = function(model, mvSaved, target, control) {
        regular_sampler <- nimbleFunctionList(sampler_BASE)
        regular_sampler[[1]] <- control$wrapped_conf$buildSampler(model, mvSaved)
        dcrpNode <- control$dcrpNode
        clusterID <- control$clusterID
    },
    run = function() {
        if(any(model[[dcrpNode]] == clusterID)) regular_sampler[[1]]$run()
    },
    methods = list(
        reset = function() {regular_sampler[[1]]$reset()}
    ))


getSamplesDPmeasureNames <- function(clusterVarInfo, model, truncG, p) {
  result <- NULL
  for(j in 1:p) {
    tildeNodesModel <- model$expandNodeNames(clusterVarInfo$clusterVars[j], returnScalarComponents=TRUE) # tilde nodes j in model
    allIndexes <- 1:length(tildeNodesModel)
    
    clusterID <- 1
    tildeNodesPerClusterID <- model$expandNodeNames(clusterVarInfo$clusterNodes[[j]][clusterVarInfo$clusterIDs[[j]] == clusterID], returnScalarComponents=TRUE) # tilde nodes in cluster with id 1
    aux <- match(tildeNodesModel, tildeNodesPerClusterID, nomatch = 0) 
    cn <-  tildeNodesModel[which(aux != 0)]
    
    tmp <- matrix('', length(cn), truncG)
    for(i in seq_along(cn)) {   # for each element of the cluster parameters, replicate 1:truncG times
      expr <- parse(text = cn[i])[[1]]
      tmp[i, ] <- sapply(seq_len(truncG),
                         function(idx) {
                           expr[[2+clusterVarInfo$indexPosition[j]]] <- as.numeric(idx)
                           return(deparse(expr))
                         })
    }
    result <- rbind(result , tmp)
  }
  return(c(result))
}
