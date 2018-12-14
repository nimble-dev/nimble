## samples from measure G after initial MCMC is run on a CRP-based model
## Used when syntax xi[1:N] ~ dCRP(conc) is used in BUGS.

##-----------------------------------------
##  Wrapper function for sampleDPmeasure
##-----------------------------------------

#' Get posterior samples for a Dirichlet process measure
#'
#' EXPERIMENTAL This function obtains posterior samples from a Dirichlet process distributed random measure of a model specified using the \code{dCRP} distribution.
#'
#' @param MCMC an MCMC class object, either compiled or uncompiled.
#' 
#' @author Claudia Wehrhahn and Christopher Paciorek
#' 
#' @export
#' @details
#' This function provides samples from a random measure having a Dirichlet process prior. Realizations are almost surely discrete and represented by a (finite) stick-breaking representation (Sethuraman, 1994), whose atoms (or point masses) are independent and identically distributed. This sampler can only be used with models containing a \code{dCRP} distribution . 
#'
#' The \code{MCMC} argument is an object of class MCMC provided by \code{buildMCMC}, or its compiled version. The MCMC should already have been run, as \code{getSamplesDPmeasure} uses the parameter samples to  generates samples for the random measure. Note that the monitors associated with that MCMC must include the cluster membership variable (which has the \code{dCRP} distribution), the cluster parameter variables, all variables directly determining the \code{dCRP} concentration parameter, and any stochastic parent variables of the cluster parameter variables. See \code{help(configureMCMC)} or \code{help(addMonitors)} for information on specifying monitors for an MCMC.
#' 
#' The truncation level of the random measure is determined based on a fixed error of approximation and by the posterior samples of the concentration parameter, if random. The error of approximation is the tail probability of the random measure (Section 4 in Ishwaran and Zarepour, 2000).
#'  
#' The returned list contains a matrix with samples from the random measure (one sample per row) and the truncation level. The stick-breaking weights are named \code{weights} and the atoms, or point masses, are named based on the cluster variables in the model.
#' 
#' @seealso \code{\link{buildMCMC}}, \code{\link{configureMCMC}}, 
#' @references
#'
#' Sethuraman, J. (1994). A constructive definition of Dirichlet priors. \emph{Statistica Sinica}, 639-650.
#'
#' Ishwaran, H., and Zarepour, M. (2000). Markov chain Monte Carlo in approximate Dirichlet and beta two-parameter process hierarchical models. \emph{Biometrika}, 87(2), 371-390.
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
getSamplesDPmeasure <- function(MCMC) {
  if(exists('model',MCMC, inherits = FALSE))
    compiled <- FALSE else compiled <- TRUE
    if(compiled) {
      if(!exists('Robject', MCMC, inherits = FALSE) || !exists('model', MCMC$Robject, inherits = FALSE))
        stop("getSamplesDPmeasure: problem with finding model object in compiled MCMC")
      model <- MCMC$Robject$model
      mvSamples <- MCMC$Robject$mvSamples
    } else {
      model <- MCMC$model
      mvSamples <- MCMC$mvSamples
    }
    rsampler <- sampleDPmeasure(model, mvSamples)
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
    
    namesAtoms <- c()
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
  
  setup=function(model, mvSaved){
    
    ## Check if the mvSaved is compiled or not.
    # mvIsCompiled <- exists('dll', envir = mvSaved)
    # if( mvIsCompiled ) {
    #  stop("sampleDPmeasure: modelValues object has to be an uncompiled object.\n")
    #}
    
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
      stop( 'sampleDPmeasure: The node having the dCRP distribution has to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    
    ## Find the cluster variables, named tildeVars
    targetElements <- model$expandNodeNames(dcrpNode, returnScalarComponents = TRUE)
    tildeVars <- nimble:::findClusterVars(model, targetElements[1]) 
    if( is.null(tildeVars) )  ## probably unnecessary as checked in CRP sampler, but best to be safe
      stop('sampleDPmeasure: The model should have at least one cluster variable.\n')
    
    ## Check that cluster variables are monitored.
    counts <- tildeVars %in% mvSavedVars
    if( sum(counts) != length(tildeVars) ) 
      stop('sampleDPmeasure: The node(s) representing the cluster variables must be monitored in the MCMC (and therefore stored in the modelValues object).\n')  
    
    ## Check that tilde nodes are continuous and univariate variables (for avoiding long trials in simulating atoms for G in run code:
    #if(any(model$isDiscrete(tildeVars)))
    #  stop('sampleDPmeasure: cluster variables should be continuous random variables.\n')
    
    #if(any(model$isMultivariate(tildeVars)))
    #  stop( 'sampleDPmeasure: only univariate cluster variables are allowed.\n' )
    
    ## Getting all stochastic parent nodes of cluster variables (needed for simulating new tildeVar values):
    parentNodesTildeVars <- NULL
    tildeVarsElements <- list()
    for(i in seq_along(tildeVars) ) {
      tildeVarsElements[[i]] <- model$expandNodeNames(tildeVars[i])
    }        
    candidateParentNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes %in% unlist(tildeVarsElements)]
    for(i in seq_along(candidateParentNodes)) {
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE)
      for(j in seq_along(tildeVars)) {
        if(sum(aux == tildeVarsElements[[j]][1]))
          parentNodesTildeVars <- c(parentNodesTildeVars, candidateParentNodes[i])
      }
    }
    if(length(parentNodesTildeVars)) {
      parentNodesTildeVarsDeps <- model$getDependencies(parentNodesTildeVars, self = FALSE)
    } else parentNodesTildeVarsDeps <- NULL
    ## make sure tilde vars are included (e.g., if a tilde var has no stoch parents) so they get simulated
    parentNodesTildeVarsDeps <- model$expandNodeNames(c(parentNodesTildeVarsDeps, tildeVars), sort = TRUE)
    
    if(!all(model$getVarNames(nodes = parentNodesTildeVars) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the cluster variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    if(is.null(parentNodesTildeVars)) parentNodesTildeVars <- tildeVars  ## to avoid NULL which causes compilation issues

      ## CLAUDIA: here is my suggestion - we first check fully IID for each tilde var and then
      ## we check the normal-invgamma special case. 
      
    # checks for iid assumption of cluster parameters in the definition of random measure G
    # i) check that all parentNodesTildeVarsDeps have same distribution  and parameters
      nonIID <- FALSE
      for(i in seq_along(tildeVars) ) {
          clusterNodes <- model$expandNodeNames(tildeVars[i])
          valueExprs <- sapply(clusterNodes, function(x) model$getValueExpr(x))
          names(valueExprs) <- NULL
          if(length(unique(valueExprs)) != 1) {
              nonIID <- TRUE
          }
      }
      if(nonIID && length(tildeVars) == 2) {  ## check for normal-invgamma conjugacy
          stochNodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
          clusterNodes1 <- model$expandNodeNames(tildeVars[1])
          clusterNodes2 <- model$expandNodeNames(tildeVars[2])
          ## Avoid non-nodes from truncated clustering,
          ## e.g., avoid 'thetaTilde[3:10]' if only have 2 thetaTilde nodes but 10 obs.
          clusterNodes1 <- clusterNodes1[clusterNodes1 %in% stochNodes]
          clusterNodes2 <- clusterNodes2[clusterNodes2 %in% stochNodes]
          exampleNodes <- c(clusterNodes1[1], clusterNodes2[1])
          dists <- c(model$getDistribution(clusterNodes1[1]), model$getDistribution(clusterNodes2[1]))
          if(dists[1] == "dinvgamma" && dists[2] ==  "dnorm") {  ## put in order so dnorm node is first
              exampleNodes <- c(clusterNodes2[1], clusterNodes1[1])
              dists <- c("dnorm", "dinvgamma")
              tmp <- clusterNodes1; clusterNodes1 <- clusterNodes2; clusterNodes2 <- tmp
          }
          if(dists[1] == "dnorm" && dists[2] == "dinvgamma") {
              ## Check for conjugacy so we know we are in dnorm_invgamma case
              conjugacy_dnorm <- model$checkConjugacy(exampleNodes[1], restrictLink = 'identity')
              conjugacy_dinvgamma <- model$checkConjugacy(exampleNodes[2])
              if(length(conjugacy_dnorm) && length(conjugacy_dinvgamma) &&
                 sum(conjugacy_dinvgamma[[1]]$control$dep_dnorm == exampleNodes[1])) {
                  conjugate <- TRUE
              }
              
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
                  if(!length(grep(clusterNodes2[i], varText, fixed = TRUE))) {
                      varExprs[[i]] <- NULL
                  } else varExprs[[i]] <- parse(text = gsub(clusterNodes2[i], "a", varText, fixed = TRUE))[[1]]
              }
              if(length(unique(varExprs)) != 1)
                  conjugate <- FALSE
              
              ## Check that cluster variance nodes are IID
              valueExprs <- sapply(clusterNodes2, function(x) model$getValueExpr(x))
              names(valueExprs) <- NULL
              if(length(unique(valueExprs)) != 1)
                  conjugate <- FALSE
          }
          if(conjugate) nonIID <- FALSE
      }

    if(nonIID) stop('sampleDPmeasure: cluster parameters have to be independent and identically distributed. \n')
    
    fixedConc <- TRUE # assume that conc parameter is fixed. This will change in the if statement if necessary
    
    ## Get all parents of xi (membership) variable (i.e., nodes involved in concentration parameter), potentially
    ## needed to determine concentration parameter for each iteration of MCMC
    #parentNodesXi <- NULL
    #candidateParentNodes <- model$getNodeNames(includeData = FALSE)
    #candidateParentNodes <- candidateParentNodes[!candidateParentNodes == dcrpNode]
    #for(i in seq_along(candidateParentNodes)){
    #  aux <- model$getDependencies(candidateParentNodes[i], self = FALSE) 
    #  if(sum(aux == dcrpNode)) 
    #    parentNodesXi <- c(parentNodesXi, candidateParentNodes[i])
    #}
    
    #if(length(parentNodesXi)) { # concentration parameter is random (or at least a function of other quantities)
    #  ## Check that values stored in mvSaved are sufficient to calculate dcrp concentration.
    #  ## Do this by creating a model containing all NAs and then copying in variables that are saved 
    #  ## and then checking getParam gives a non-NA.
    #  fixedConc <- FALSE
      
    #  ## Which parent nodes are saved.
    #  parentNodesXi <- parentNodesXi[parentNodesXi %in% mvSavedVars]  
      
    #  ## Create model with NA values
    #  verbosity <- nimbleOptions('verbose')
    #  nimbleOptions(verbose = FALSE)
    #  modelWithNAs <- model$modelDef$newModel(check = FALSE, calculate = FALSE)
    #  nimbleOptions(verbose = verbosity)
      
      ## Note that this check could fail if current state of model has NAs that result in NA in conc parameter,
      ## but we don't want to use mvSaved as MCMC may not have been run at this point.
      ## Need to think more about this.
    #  nimCopy(from = model, to = modelWithNAs, nodes = parentNodesXi) 
    #  if(length(parentNodesXi)) {
    #    ## These nodes need to be calculated to determine conc param in run code
    #    parentNodesXiDeps <- model$getDependencies(parentNodesXi, self = FALSE, determOnly = TRUE)  ## excludes dcrpNode
    ##    modelWithNAs$calculate(parentNodesXiDeps)
    #    dcrpParam <- modelWithNAs$getParam(dcrpNode, 'conc')
    #    if(is.na(dcrpParam)) 
    #      stop('sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n') 
    #  } else stop( 'sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n') 
    #} else { ## placeholder since parentNodesXi must only have nodes in the mvSaved for correct compilation
    #  parentNodesXiDeps <- dcrpNode
    #  parentNodesXi <- dcrpNode # also used in run code
    #}  
    
    
    # test this version of doing things....
    parentNodesXi <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE, stochOnly = TRUE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes == dcrpNode]
    for(i in seq_along(candidateParentNodes)) {
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE)
      if(sum(aux == dcrpNode)) {
        parentNodesXi <- c(parentNodesXi, candidateParentNodes[i])
      }
    }
    
    if(length(parentNodesXi)) {
      fixedConc <- FALSE
      parentNodesXiDeps <- model$getDependencies(parentNodesXi, self = FALSE)
      parentNodesXiDeps <- parentNodesXiDeps[!parentNodesXiDeps == dcrpNode]
    } else {
      parentNodesXiDeps <- dcrpNode
    }
    
    if(!all(model$getVarNames(nodes = parentNodesXi) %in% mvSavedVars))
      stop('sampleDPmeasure: The stochastic parent nodes of the membership variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    if(is.null(parentNodesXi)) parentNodesXi <- dcrpNode  ## to avoid NULL which causes compilation issues
    
    
    
    
    dataNodes <- model$getDependencies(targetElements[1], stochOnly = TRUE, self = FALSE)
    N <- length(model$getDependencies(targetElements, stochOnly = TRUE, self = FALSE))
    p <- length(tildeVars)
    lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE))
    nTilde <- numeric(p+1)
    dimTildeNim <- numeric(p+1) # nimble dimension (0 is scalar, 1 is 2D array, 2 is 3D array)
    dimTilde <- numeric(p+1) # dimension to be used in run code
    #dimTildeNimAux <- numeric(p) #
    for(i in 1:p) {
      elementsTildeVars <- model$expandNodeNames(tildeVars[i], returnScalarComponents = TRUE)
      dimTildeNim[i] <- model$getDimension(elementsTildeVars[1]) # if the node is deterministic returns NA. I should get the dimension of the parent node?
      dimTilde[i] <- lengthData^(dimTildeNim[i]) 
      nTilde[i] <- length(values(model, tildeVars[i])) / dimTilde[i]
    }
    if(any(nTilde[1:p] != nTilde[1])){
      stop('sampleDPmeasure: All cluster parameters must have the same number of parameters.\n')
    }
    
    
    ## The error of approximation G is given by (conc / (conc +1))^{truncG-1}. 
    ## we are going to define an error of aproximation and based on the posterior values of the conc parameter define the truncation level of G
    ## the error is between errors that are considered very very small in the folowing papers
    ## Ishwaran, H., & James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. Journal of the American Statistical Association, 96(453), 161-173.
    ## Ishwaran, H., & Zarepour, M. (2000). Markov chain Monte Carlo in approximate Dirichlet and beta two-parameter process hierarchical models. Biometrika, 87(2), 371-390.
    approxError <- 1e-15
    
    ## Storage object to be sized in run code based on MCMC output (Claudia note change to comment)
    samples <- matrix(0, nrow = 1, ncol = 1)   
    ## Tuncation level of the random measure 
    truncG <- 0 
    
    setupOutputs(lengthData)
  },
  
  run=function(){
    
    niter <- getsize(mvSaved) # number of iterations in the MCMC
    
    # defining the truncation level of the random measure's representation:
    if( fixedConc ) {
      dcrpAux <- model$getParam(dcrpNode, 'conc')
      concSamples <- nimNumeric(length = niter, value = dcrpAux)
    } else {
      concSamples <- numeric(niter)
      for( iiter in 1:niter ) {
        nimCopy(from = mvSaved, to = model, nodes = parentNodesXi, row=iiter) 
        model$calculate(parentNodesXiDeps)
        concSamples[iiter] <- model$getParam(dcrpNode, 'conc')
      }
      dcrpAux <- mean(concSamples)
    }
    
    truncG <<- log(approxError) / log(dcrpAux / (dcrpAux+1)) + 1
    truncG <<- round(truncG)
    #approxError <- (dcrpAux / (dcrpAux +1))^(truncG-1)
    # I think is good to send message indicating what the truncation level is for an approximation error smaller than to 10^(-10)
    # nimCat('sampleDPmeasure: Approximating the random measure by a finite stick-breaking representation with an error smaller than 1e-10, leads to a truncation level of ', truncG, '.\n')
    
    ## Storage object: matrix with nrow = number of MCMC iterations, and ncol = (1 + p)*truncG, where
    ## truncG the truncation level of the random measure G (an integer given by the values of conc parameter)
    ## (p+1) denoted the number of parameters, each of length truncG, to be stored: p is the number of cluster components (length(tildeVars)) and 1 is for the weights.
    samples <<- matrix(0, nrow = niter, ncol = truncG*(sum(dimTilde)+1)) 
    
    for(iiter in 1:niter){
      checkInterrupt()
      
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
              if(dimTildeNim[j] == 0) { # scalars
                uniqueValues[index, jcol] <- values(model, tildeVars[j])[dimTilde[j]*(range[i] - 1) + l]  #   
              }
              if(dimTildeNim[j] == 2) { #  matrices
                uniqueValues[index, jcol] <- values(model, tildeVars[j])[dimTilde[j]*(range[i] - 1) + l]  #   
              }
              if(dimTildeNim[j] == 1) { # vectors
                uniqueValues[index, jcol] <- values(model, tildeVars[j])[range[i] + (l-1)*nTilde[j]]  #   
              }
              jcol <- jcol + 1
            }
          }
          index <- index+1
        }
      }
      probs[index] <- concSamples[iiter] 
      newValueIndex <- index 
      
      ## computing G(.) = sum_{l=1}^{truncG} w_l delta_{atom_l} (.):
      vaux <- rbeta(1, 1, concSamples[iiter] + N)
      v1prod <- 1
      Taux <- 1
      paramAux <- numeric(sum(dimTilde))
      
      ## copy tilde parents into model for use in simulation below when simulate from G_0  
      nimCopy(mvSaved, model, parentNodesTildeVars, row = iiter)
      
      ## first sampled values: w_1 and atom_1
      index <- rcat(prob = probs[1:newValueIndex])
      if(index == newValueIndex){   # sample from G_0
        model$simulate(parentNodesTildeVarsDeps)
        sumCol <- truncG
        for(j in 1:p){
          for(l in 1:dimTilde[j]) {
            if(dimTildeNim[j] == 0) {
              samples[iiter, sumCol + dimTilde[j]*(Taux - 1) +l] <<- values(model, tildeVars[j])[l]  #   
            }
            if(dimTildeNim[j] == 2) {
              samples[iiter, sumCol + dimTilde[j]*(Taux - 1) +l] <<- values(model, tildeVars[j])[l]  #   
            }
            if(dimTildeNim[j] == 1) {
              samples[iiter, sumCol + dimTilde[j]*(Taux - 1) +l] <<- values(model, tildeVars[j])[(l-1)*nTilde[j] + 1]  # 
            }
            
          }
          sumCol <- sumCol + dimTilde[j]*truncG 
        }
      } else {   # sample one of the existing values
        sumCol <- truncG
        jcol <- 1
        for(j in 1:p){
          for(l in 1:dimTilde[j]) {
            samples[iiter, sumCol + dimTilde[j]*(Taux - 1) +l] <<- uniqueValues[index, jcol]   # 
            jcol <- jcol + 1
          }
          sumCol <- sumCol + dimTilde[j]*truncG 
        }
      }
      samples[iiter, Taux] <<- vaux 
      Taux <- Taux + 1
      
      # the rest of the values: w_l and atom_l, l=1, ..truncG-1
      while(Taux <= truncG){
        index <- rcat(prob = probs[1:newValueIndex])
        if(index == newValueIndex){  # sample from G_0
          model$simulate(parentNodesTildeVarsDeps)
          jcol <- 1
          for(j in 1:p){
            for(l in 1:dimTilde[j]) {
              if(dimTildeNim[j] == 0 ) {
                paramAux[jcol] <- values(model, tildeVars[j])[l]   #   
              }
              if(dimTildeNim[j] == 2 ) {
                paramAux[jcol] <- values(model, tildeVars[j])[l]   #   
              }
              if(dimTildeNim[j] == 1) {
                paramAux[jcol] <- values(model, tildeVars[j])[(l-1)*nTilde[j] +1] 
              }
              jcol <- jcol + 1
            }
          }
        } else{  # sample one of the existing values
          jcol <- 1
          for(j in 1:p){
            for(l in 1:dimTilde[j]) {
              paramAux[jcol] <- uniqueValues[index, jcol]   # 
              jcol <- jcol + 1
            }
          }
        }
        condaux <- samples[iiter, truncG + seq(1, dimTilde[1]*(Taux - 1), by=dimTilde[1]) ] == paramAux[1]  # 1:(dimTilde[1]*(Taux - 1)) check if we sample a new atom or an atom that is in G already
        if(sum(condaux) > 0) { # the atom already exists and we have to update the weights and not include a new value of the params 
          repindex = 1
          while(!condaux[repindex]){
            repindex = repindex + 1
          }
          v1prod <- v1prod * (1-vaux)
          vaux <- rbeta(1, 1, concSamples[iiter]+N)
          samples[iiter, repindex] <<- samples[iiter, repindex] + vaux * v1prod
        } else { # augment the truncation and keep the same parameters
          sumCol <- truncG
          jcol <- 1
          for(j in 1:p){
            for(l in 1:dimTilde[j]) {
              samples[iiter, sumCol + dimTilde[j]*(Taux - 1) +l] <<- paramAux[jcol]   # 
              jcol <- jcol + 1
            }
            sumCol <- sumCol + dimTilde[j]*truncG 
          }
          v1prod <- v1prod * (1-vaux)
          vaux <- rbeta(1, 1, concSamples[iiter]+N)
          if(Taux != truncG) {
            samples[iiter, Taux] <<- vaux * v1prod   
          } else {
            samples[iiter, Taux] <<- 1 * v1prod 
          }
          Taux <- Taux + 1
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
    N <- length(model[[xiNode]])
  },
  
  
  run = function() {
    shapeParam <- model$getParam(target, 'shape')
    rateParam <- model$getParam(target, 'rate')
    conc <- model[[target]] 
    xi <- model[[xiNode]]
    
    occupied <- numeric( N)
    for( i in 1:N )
      occupied[xi[i]] <- 1
    k <- sum(occupied)
    
    aux1 <- shapeParam + k
    aux2 <- aux1 - 1
    
    ## generating augmented r.v. and computing the weight.
    x <- rbeta(1, conc+1, N)
    aux3 <- rateParam - log(x)
    w <- aux2/(aux2 + N*aux3)
    
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



#-----------------------------------
# Conjugate cases in Sampler CRP
#-----------------------------------
## we need a base class because it is possible (perhaps unlikely) that
## a user model might have two uses of dCRP samplers that are different versions
## e.g., a nonconjugate and a dnorm_dnorm conjugate
## to allow for this we need to use a one-element nimbleFunctionList that uses
## this virtual base class
CRP_helper <- nimbleFunctionVirtual(
  methods = list(
    storeParams = function() { },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
    },
    sample = function(i = integer(), j = integer()) {}
  )
)

CRP_nonconjugate <- nimbleFunction(
  name = "CRP_nonconjugate",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
  },
  methods = list(
    storeParams = function() {},  ## nothing needed for non-conjugate
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      return(model$getLogProb(dataNodes[i]))
    },
    sample = function(i = integer(), j = integer() ) { ## sample from the base measure of the DP
      p <- length(marginalizedVar)
      nTildeVars <- length(marginalizedNodes) / p
      for(l in seq_along(marginalizedVar)) {
        model$simulate(marginalizedNodes[(l-1)*nTildeVars + j])
      }
    }
  )
)

CRP_conjugate_dnorm_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dnorm_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      #model[[marginalizedVar]][j] <<- rnorm(1, postMean, sqrt(postVar)) 
      values(model, marginalizedNodes[j]) <<- rnorm(1, postMean, sqrt(postVar)) 
    }
  )
)


CRP_conjugate_dnorm_invgamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dnorm_invgamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
    priorMean <- nimNumeric(1)
    kappa0 <- nimNumeric(1)
    priorShape <- nimNumeric(1)
    priorScale <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priorMean <<- model$getParam(marginalizedNodes[1], 'mean_location')
      kappa0 <<- model$getParam(marginalizedNodes[1], 'mean_scale')
      priorShape <<- model$getParam(marginalizedNodes[1], 'var_shape')
      priorScale <<- model$getParam(marginalizedNodes[1], 'var_scale')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      c1 <- priorShape * log(priorScale) + lgamma(priorShape + 1/2) + log(kappa0) -
       lgamma(priorShape) - log(2) - log(pi) - log(1 + kappa0)
      c2 <- - (priorShape  + 1/2) * (priorScale + kappa0 * (y - priorMean)^2 / (2*(1+kappa0)) )
      return(c1 + c2)
    },
    sample = function(i = integer(), j = integer()) {
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][2] <<- rinvgamma(shape = priorShape + 1/2,
                                      scale = priorScale + kappa0 * (y - priorMean)^2 / (2*(1+kappa0)) )
      model[[marginalizedVar]][1] <<- rnorm(1, (kappa0 * priorMean + y)/(1 + kappa0), 
                                            sd = sqrt(model[[marginalizedVar]][2] / (1+kappa0))) 
    }
  )
)


CRP_conjugate_dgamma_dpois <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dpois",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      values(model, marginalizedNodes[j]) <<- rgamma(1, shape = priorShape + y, rate = priorRate + 1)
    }
  )
)


CRP_conjugate_dgamma_dnorm <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dnorm",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      return(-0.5*log(2*pi) + priorShape * log(priorRate) - lgamma(priorShape) -
               lgamma(priorShape + 0.5) - (priorShape + 0.5)*log(priorRate + (y-dataMean)^2/2))
    },
    sample = function(i = integer(), j = integer()) {
      dataMean <- model$getParam(dataNodes[i], 'mean')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- rgamma(1, shape = priorShape + 0.5, rate = priorRate + (y-dataMean)^2/2)
    }
  )
)



CRP_conjugate_dbeta_dbern <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dbern",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      values(model, marginalizedNodes[j]) <<- rbeta(1, shape1=priorShape1+y, shape2=priorShape2+1-y)
    }
  )
)

CRP_conjugate_dbeta_dbin <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dbin",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      values(model, marginalizedNodes[j]) <<- rbeta(1, shape1=priorShape1+y, shape2=priorShape2+dataSize-y)
    }
  )
)



CRP_conjugate_dbeta_dnegbin <- nimbleFunction(
  name = "CRP_conjugate_dbeta_dnegbin",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      values(model, marginalizedNodes[j]) <<- rbeta(1, shape1=priorShape1+dataSize, shape2=priorShape2+y)
    }
  )
)

CRP_conjugate_dgamma_dexp <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dexp",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      values(model, marginalizedNodes[j]) <<- rgamma(1, shape=priorShape+1, rate=priorRate+y)
    }
  )
)


CRP_conjugate_dgamma_dgamma <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dgamma",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      values(model, marginalizedNodes[j]) <<- rgamma(1, shape=datashape+priorShape, rate=priorRate+y)
    }
  )
)


CRP_conjugate_dgamma_dweib <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dweib",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
      values(model, marginalizedNodes[j]) <<- rgamma(1, shape=1+priorShape, rate=priorRate+y^dataShape)
    }
  )
)


CRP_conjugate_dgamma_dinvgamma <- nimbleFunction(
  name = "CRP_conjugate_dgamma_dinvgamma",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
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
                lgamma(dataShape) - lgamma(dataShape) - (dataShape + priorShape)*log(priorRate + 1/y))
    },
    sample = function(i = integer(), j = integer()) {
      dataShape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      values(model, marginalizedNodes[j]) <<- rgamma(1, shape=dataShape+priorShape, rate=priorRate+1/y)
    }
  )
)


CRP_conjugate_ddirch_dmulti <- nimbleFunction(
  name = "CRP_conjugate_ddirch_dmulti",
  contains = CRP_helper,
  setup = function(model, marginalizedVar, marginalizedNodes, dataNodes) {
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
      return(lfactorial(n) - sum(lfactorial(y)+lgamma(priorAlpha))+
               lgamma(sum(priorAlpha)) + sum(lgamma(priorAlpha+y)) - lgamma(sum(priorAlpha+y)))
    },
    sample = function(i = integer(), j = integer()) {
      y <- values(model, dataNodes[i])
      model[[marginalizedVar]][j, ] <<- rdirch(alpha=priorAlpha+y)
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
    ## note that even in inefficient case, we need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.
    
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents = TRUE)
    targetVar <- model$getVarNames(nodes = target)  
    n <- length(targetElements) 
    
    # first check that the sampler can be used: we need one observation per random index
    nObs <- length(model$getDependencies(targetElements, stochOnly = TRUE, self = FALSE))
    if(n != nObs)
      stop("sampler_CRP: The length of membership variable and observations has to be the same.\n")

    clusterVarInfo <- findClusterVars(model, targetElements[1], returnIndexInfo = TRUE)
    tildeVars <- clusterVarInfo$clusterVars
    if(is.null(tildeVars))
        stop('sampler_CRP:  The model should have at least one cluster variable.\n')
    ## Check for indexing that would cause errors in using sample() with 'j' as the index of the set of cluster nodes, e.g., mu[xi[i]+2] or mu[some_function(xi[i])]
    ## CLAUDIA:  
    if(any(clusterVarInfo$targetInFunction))
        stop("sampler_CRP: At the moment, NIMBLE's CRP MCMC sampling requires simple indexing and cannot work with indexing such as 'mu[xi[i]+2]'.")

    calcNodes <- unique(c(calcNodes, model$getDependencies(tildeVars)))

    # check that the number of parameters in cluster parameters (if more than one) are the same  (multivariate case considered)
    dataNodes <- model$getDependencies(targetElements[1], stochOnly = TRUE, self = FALSE) 
    lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE)) # for vector data gives its length 
    p <- length(tildeVars)
    nTilde <- numeric(p)
    for(i in 1:p) {
      elementsTildeVars <- model$expandNodeNames(tildeVars[i], returnScalarComponents = TRUE)
      dimTildeNim <- model$getDimension(elementsTildeVars[i])
      nTilde[i] <- length(values(model, tildeVars[i])) / (lengthData)^dimTildeNim
    }
    # in the univariate case: nTilde <- sapply(tildeVars, function(x) length(model[[x]]))
    if(any(nTilde != nTilde[1])){
      stop('sampleDPmeasure: All cluster parameters must have the same number of observations.\n')
    }
    
    if(length(unique(nTilde)) != 1)
      stop('sampler_CRP: In a model with multiple cluster parameters, the number of those parameters must all be the same.\n')
    
    min_nTilde <- min(nTilde) ## we need a scalar for use in run code
    if(min_nTilde < n)
      warning('sampler_CRP: The number of cluster parameters is less than the number of potential clusters. The MCMC is not strictly valid if it ever proposes more components than cluster parameters exist; NIMBLE will warn you if this occurs.\n')
    
    
    # determine if concentration parameter is fixed or random (code similar to the one in sampleDPmeasure function):
    fixedConc <- TRUE
    stochNodes <- model$getNodeNames(stochOnly = TRUE)
    distributions <- model$getDistribution(stochNodes) 
    dcrpIndex <- which(distributions == 'dCRP')
    dcrpNode <- stochNodes[dcrpIndex] 
    
    parentNodesXi <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes == dcrpNode]
    for(i in seq_along(candidateParentNodes)){
      tmp <- model$getDependencies(candidateParentNodes[i], self = FALSE) 
      if(sum(tmp == dcrpNode)) 
        parentNodesXi <- c(parentNodesXi, candidateParentNodes[i])
    }
    if(length(parentNodesXi)) {
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
        stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self = FALSE) 
        detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
        if(length(stochDeps) != 1) 
          stop("sampler_CRP: NIMBLE cannot currently assign a sampler to a dCRP node unless each membership element is associated with a single observation.\n")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
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

    
    helperFunctions <- nimble:::nimbleFunctionList(CRP_helper)
    
    ## use conjugacy to determine which helper functions to use
    conjugacyResult <- checkCRPconjugacy(model, target)
    if(is.null(conjugacyResult)) {
      sampler <- 'CRP_nonconjugate'
    } else 
      sampler <- switch(conjugacyResult,
                        conjugate_dnorm_dnorm = 'CRP_conjugate_dnorm_dnorm',
                        #conjugate_dnorm_invgamma_dnorm = 'CRP_conjugate_dnorm_invgamma_dnorm',
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

    tildeNodes <- model$expandNodeNames(tildeVars, sort = TRUE)
    ## check that tildeVars are grouped together by variable because sample() assumes this.
    ## E.g., c('sigma','sigma','mu','mu') for n=2,p=2 is ok but c('sigma','mu','sigma','mu') is not.
    ## CLAUDIA:  
    if(p > 1) {  
        tmp <- matrix(tildeNodes, ncol = p)
        if(any(sapply(tmp, function(x) length(unique(x)) > 1)))
            stop("sampler_CRP: Current MCMC CRP sampling cannot handle the structure of the multiple cluster variables.")
    }
      
    helperFunctions[[1]] <- eval(as.name(sampler))(model, tildeVars, tildeNodes, dataNodes)
    
    curLogProb <- numeric(n)
  },
  
  
  run = function() {
    
    isNonParam <- TRUE # indicates if the MCMC sampling is nonparametric or not. Is nonparametric when there are more cluster parameters than required clusters.
    conc <- model$getParam(target, 'conc')
    helperFunctions[[1]]$storeParams()
    
    xi <- model[[target]]
    
    ## Find unique values in model[[target]]. I'm not relabeling the unique values.
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
    while(xiCounts[kNew] > 0 & kNew < n) { # need to make sure don't go beyond n
      kNew <- kNew + 1
    }
    if( kNew == n & xiCounts[kNew] > 0 ) { # case  xi=1:n
      kNew <- 0
    }
    if(kNew > min_nTilde & min_nTilde != n) {
      if(fixedConc) {
        nimCat('CRP_sampler: This MCMC is for a parametric model. The MCMC attempted to use more components than the number of cluster parameters. To have a sampler for a nonparametric model increase the number of cluster parameters.\n')
      } else {
        nimCat('CRP_sampler: This MCMC is not for a proper model. The MCMC attempted to use more components than the number of cluster parameters. Please increase the number of cluster parameters.\n')
      }
      kNew <- 0 
      isNonParam <- FALSE
    }
    
    
    for(i in 1:n) { # updates one cluster membership at the time , i=1,...,n
      
      xi <- model[[target]]
      xiCounts[xi[i]] <- xiCounts[xi[i]] - 1
      
      # computing sampling probabilities and sampling an index:
      if( xiCounts[xi[i]] == 0 ) { # cluster is a singleton. First, compute probability of sampling an existing label.
        # Second, compute probability of sampling a new cluster, here, new cluster is the current cluster!
        reorderXiUniques <- numeric(min_nTilde) # here we save reordered version of xiUniques when there is a singleton. This is used latter for updating xiUniques if a component is deleted
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

        # sample and index from 1 to k
        index <- rcat( n=1, exp(curLogProb[1:k]-max(curLogProb[1:k])) )
        if(index == k) {
          newLab <- xi[i] 
          newLabCond <- TRUE
        } else {
          newLab <- reorderXiUniques[index]
          newLabCond <- FALSE
        }
        
      } else { # cluster is not a singleton. First, compute probability of sampling an existing label.
        # Second, compute probability of sampling a new cluster dependening on the value of kNew
        for(j in 1:k) { # probability of sampling an existing label
          model[[target]][i] <<- xiUniques[j] 
          if(type == 'indivCalcs') {
            if(nInterm >= 1) model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
          curLogProb[j] <<- log(xiCounts[xiUniques[j]]) + model$getLogProb(dataNodes[i]) 
        }
        
        if(kNew == 0) { # no new cluster can be created 
          curLogProb[k+1] <<- log(0)
        } else { # a new cluster can be created
          model[[target]][i] <<- kNew
          if(sampler == 'CRP_nonconjugate'){
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
      
      # updating objects
      model[[target]][i] <<- newLab
      
      if( newLabCond ) { # a component is created. It can really create a new component or keep the current label if xi_i is a singleton
        if(isNonParam & sampler != 'CRP_nonconjugate') { # updating the cluster parameters of the new cluster
          helperFunctions[[1]]$sample(i, model[[target]][i])
        }
        if( xiCounts[xi[i]] != 0) { # a component is really created
          k <- k + 1
          xiUniques[k] <- newLab 
          # updating kNew:
          kNew <- kNew + 1
          mySum <- sum(xi == kNew) 
          while(mySum > 0 & kNew < (n+1)) { # need to make sure don't go beyond length of vector
            kNew <- kNew+1
            mySum <- sum(xi == kNew)
          }
          if(kNew > min_nTilde) {
            if(fixedConc) {
              nimCat('CRP_sampler: This MCMC is for a parametric model. The MCMC attempted to use more components than the number of cluster parameters. To have a sampler for a nonparametric model increase the number of cluster parameters.\n')
            } else {
              nimCat('CRP_sampler: This MCMC is not for a proper model. The MCMC attempted to use more components than the number of cluster parameters. Please increase the number of cluster parameters.\n')
            }
            kNew <- 0
            isNonParam <- FALSE
          }
        }
        xiCounts[model[[target]][i]] <- 1
      } else { # an existing label is sampled
        if( xiCounts[xi[i]] == 0 ) { # xi_i is a singleton, a component was deleted
          k <- k - 1
          xiUniques <- reorderXiUniques
          if( kNew == 0 ) { # the sampler was not nonparametric or xi=1:n
            kNew <- xi[i] # 
            isNonParam <- TRUE # now the sampler is nonparametric if it was not
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
    }
  ), where = getLoadingNamespace()
)




#' @rdname samplers
#' @export
sampler_CRP_old <- nimbleFunction(
  name = 'sampler_CRP_old',
  contains=sampler_BASE,
    
  setup=function(model, mvSaved, target, control){
    ## note that even in inefficient case, we need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.
    
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
    conjugacyResult <- checkCRPconjugacy(model, target)
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


findClusterVars <- function(model, target, returnIndexInfo = FALSE) {
    targetVar <- model$getVarNames(nodes = target)
    clusterVars <- NULL
    if(returnIndexInfo) {
        ## This information is useful in checkCRPconjugacy in determining if we can handle conjugacy
        ## under various complicated indexing situations, such as 'p[2:4, xi[i]]'.
        numIndexes <- NULL
        indexPosition <- NULL
        indexExpr <- list()
        targetInFunction <- NULL
    }
    idx <- 1
    deps <- model$getDependencies(target, self = FALSE)
    for(dep in deps) { 
        expr <- cc_getNodesInExpr(model$getValueExpr(dep)) 
        for(j in seq_along(expr)) {
            ## look for cases like thetatilde[xi[i]] to identify 'xi' and extract 'thetaTilde'
            subExpr <- parse(text = expr[j])[[1]]
            len <- length(subExpr)
            if(len >= 3 && is.call(subExpr) && subExpr[[1]] == '[') { 
                ##  Find the cluster variables, named clusterVars, in presence of univariate and multivariate cluster parameters:
                k <- 3
                foundTarget <- FALSE
                while(k <= len && !foundTarget) {
                    if(sum(all.vars(subExpr[[k]]) == targetVar))
                        foundTarget <- TRUE else k <- k + 1
                }  
                if(foundTarget) {
                    clusterVars[idx] <- deparse(subExpr[[2]])
                    if(returnIndexInfo) {
                        indexPosition[idx] <- k-2
                        numIndexes[idx] <- len - 2
                        indexExpr[[idx]] <- subExpr
                        targetInFunction[idx] <- length(subExpr[[k]]) >= 3 &&
                            (length(subExpr[[k]]) > 3 || subExpr[[k]][[1]] != '[' ||
                             subExpr[[k]][[2]] != targetVar)
                    }
                    idx <- idx+1 
                }
            }
        }
    }
    if(returnIndexInfo) {
        return(list(clusterVars = clusterVars, numIndexes = numIndexes,
                    indexPosition = indexPosition, indexExpr = indexExpr,
                    targetInFunction = targetInFunction))
        } else return(clusterVars)
}
    
