## samples from measure G after initial MCMC is run on a CRP-based model
## Used when syntax xi[1:N] ~ dCRP(conc) is used in BUGS.

##-----------------------------------------
##  Wrapper function for sampleDPmeasure
##-----------------------------------------

#' Get posterior samples for a Dirichlet process measure
#'
#' EXPERIMENTAL This function obtains samples from the estimated Dirichlet process measure for models specified using the \code{dCRP} distribution.
#'
#' @param MCMC an MCMC class object, either compiled or uncompiled.
#' 
#' @author Claudia Wehrhahn and Christopher Paciorek
#' 
#' @export
#' @details
#'  This function provides samples from a truncated approximation to the random measure associated with the mixing distribution of a Dirichlet process mixture model. The random measure is represented by a stick-breaking representation (Sethuraman, 1994). This sampler can only be used with models containing a \code{dCRP} distribution. 
#'
#' The \code{MCMC} argument is an object of class MCMC provided by \code{buildMCMC}, or its compiled version. The MCMC should already have been run, as \code{getSamplesDPmeasure} uses the parameter samples to  generates samples for the random measure. Note that the monitors associated with that MCMC must include the cluster membership variable (which has the \code{dCRP} distribution), the cluster parameter variables, all variables directly determining the \code{dCRP} concentration parameter, and any stochastic parent variables of the cluster parameter variables. See \code{help(configureMCMC)} or \code{help(addMonitors)} for information on specifying monitors for an MCMC.
#' 
#' The truncation level of the random measure is determined based on a fixed error of approximation and by the posterior samples of the concentration parameter, if random. The error of approximation is the tail probability of the random measure (Section 4 in Ishwaran and Zarepour, 2000).
#'  
#' The returned object is a matrix containing samples from the truncated approximation of the random measure (one row per sample), with columns for the weights and the cluster variables. The stick-breaking weights are named \code{weights} and the atoms, or point masses, are named based on the cluster variables in the model.
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
#'   Rmcmc <- buildMCMC(conf)
#'   Cmodel <- compileNimble(model)
#'   Cmcmc <- compileNimble(Rmcmc, project = model)
#'   runMCMC(Cmcmc, niter = 1000)
#'   outputG <- getSamplesDPmeasure(Cmcmc)
#'   samples <- outputG$samples
#'   truncation <- output$trunc
#' }
getSamplesDPmeasure <- function(MCMC) {
  if(exists('model', MCMC))
    compiled <- FALSE else compiled <- TRUE
    if(compiled) {
      if(!exists('Robject', MCMC) || !exists('model', MCMC$Robject))
        stop("getSamplesDPmeasure: problem with finding model object in compiled MCMC")
      model <- MCMC$Robject$model
      dataNodes <- model$getNodeNames(dataOnly = TRUE) 
      lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE))
      mvSamples <- MCMC$Robject$mvSamples
    } else {
      model <- MCMC$model
      dataNodes <- model$getNodeNames(dataOnly = TRUE) 
      lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE))
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
        for(k in 1:lengthData){
          for(l in 1:truncG) {
            namesAtoms[inames] <- paste0(namesVars[j], "[", l, ",", k, "]")
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
    mvIsCompiled <- exists('dll', envir = mvSaved)
    if( mvIsCompiled ) {
      stop("sampleDPmeasure: modelValues object has to be an uncompiled object.\n")
    }
    
    ## Determine variables in the mv object and nodes/variables in the model.
    mvSavedVars <- mvSaved$varNames
    
    stochNodes <- model$getNodeNames(stochOnly = TRUE)
    dataNodes <- model$getNodeNames(dataOnly = TRUE) 
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
    tildeVars <- NULL
    itildeVar <- 1
    dep <- model$getDependencies(targetElements[1], self = FALSE)
    for(i in seq_along(dep)) { 
      expr <- cc_getNodesInExpr(model$getValueExpr(dep[i])) 
      for(j in seq_along(expr)) {
        ## look for cases like thetatilde[xi[i]] to identify 'xi' and extract 'thetaTilde'
        tmpexpr <- parse(text = expr[j])[[1]]
        ltmpexpr <- length(tmpexpr)
        if(ltmpexpr >= 3 && is.call(tmpexpr) && tmpexpr[[1]] == '[') { 
          #  Find the cluster variables, named tildeVars, in presence of univariate and multivariate cluster parameters:
          k <- 1
          foundTarget <- FALSE
          while(k <= ltmpexpr && foundTarget == FALSE) {
            foundTarget <- all.vars(tmpexpr[[k]]) == dcrpVar
            if(sum(foundTarget) == 0) {  # to avoid having foundTarget equal to logical(0). what other type could foundTarget be?
              foundTarget <- FALSE
            }
            k <- k + 1
          }  
          if( length(foundTarget) > 0 && sum(foundTarget) > 0 ) {
            tildeVars[itildeVar] <- deparse(tmpexpr[[2]])
            itildeVar <- itildeVar+1 
          }
        }
      }
    }
    
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
      # ToDo: add sanity check that atom are iid sampled
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
    
    fixedConc <- TRUE # assume that conc parameter is fixed. This will change in the if statement if necessary
    
    ## Get all parents of xi (membership) variable (i.e., nodes involved in concentration parameter), potentially
    ## needed to determine concentration parameter for each iteration of MCMC
    parentNodesXi <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes == dcrpNode]
    for(i in seq_along(candidateParentNodes)){
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE) 
      if(sum(aux == dcrpNode)) 
        parentNodesXi <- c(parentNodesXi, candidateParentNodes[i])
    }
    
    if(length(parentNodesXi)) { # concentration parameter is random (or at least a function of other quantities)
      ## Check that values stored in mvSaved are sufficient to calculate dcrp concentration.
      ## Do this by creating a model containing all NAs and then copying in variables that are saved 
      ## and then checking getParam gives a non-NA.
      fixedConc <- FALSE
      
      ## Which parent nodes are saved.
      parentNodesXi <- parentNodesXi[parentNodesXi %in% mvSavedVars]  
      
      ## Create model with NA values
      verbosity <- nimbleOptions('verbose')
      nimbleOptions(verbose = FALSE)
      modelWithNAs <- model$modelDef$newModel(check = FALSE, calculate = FALSE)
      nimbleOptions(verbose = verbosity)
      
      ## Note that this check could fail if current state of model has NAs that result in NA in conc parameter,
      ## but we don't want to use mvSaved as MCMC may not have been run at this point.
      ## Need to think more about this.
      nimCopy(from = model, to = modelWithNAs, nodes = parentNodesXi) 
      if(length(parentNodesXi)) {
        ## These nodes need to be calculated to determine conc param in run code
        parentNodesXiDeps <- model$getDependencies(parentNodesXi, self = FALSE, determOnly = TRUE)  ## excludes dcrpNode
        modelWithNAs$calculate(parentNodesXiDeps)
        dcrpParam <- modelWithNAs$getParam(dcrpNode, 'conc')
        if(is.na(dcrpParam)) 
          stop('sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n') 
      } else stop( 'sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n') 
    } else { ## placeholder since parentNodesXi must only have nodes in the mvSaved for correct compilation
      parentNodesXiDeps <- dcrpNode
      parentNodesXi <- dcrpNode # also used in run code
    }  
    
    N <- length(dataNodes)
    p <- length(tildeVars)
    lengthData <- length(model$expandNodeNames(dataNodes[1], returnScalarComponents = TRUE))
    nTilde <- c()
    dimTildeNim <- c() # nimble dimension (0 is scalar, 1 is 2D array, 2 is 3D array)
    dimTilde <- c() # dimension to be used in run code
    for(i in 1:p) {
      elementsTildeVars <- model$expandNodeNames(tildeVars[i], returnScalarComponents = TRUE)
      dimTildeNim[i] <- model$getDimension(elementsTildeVars[i])[[1]]
      dimTilde[i] <- lengthData^(model$getDimension(elementsTildeVars[i])[[1]]) 
      nTilde[i] <- length(values(model, tildeVars[i])) / (lengthData)^dimTildeNim[i]
    }
    if(any(nTilde != nTilde[1])){
      stop('sampleDPmeasure: All cluster parameters must have the same number of observations.\n')
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
  },
  
  run=function(){
    
    niter <- getsize(mvSaved) # number of iterations in the MCMC
    for(i in 1:p) { # so dimTildeNim can be obtained in the wrapper function
      dimTildeNim[i] <<- dimTildeNim[i]
    }
    
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
      dcrpAux <-  mean(concSamples) # quantile(concSamples, 0.9)
    }
    
    truncG <<- log(approxError) / log(dcrpAux / (dcrpAux+1)) + 1
    truncG <<- round(truncG)
    #approxError <- (dcrpAux / (dcrpAux +1))^(truncG-1)
    # I think is good to send message indicating what the truncation level is for an approximation error smaller than to 10^(-10)
    nimCat('sampleDPmeasure: Approximating the random measure by a finite stick-breaking representation with an error smaller than 1e-10, leads to a truncation level of ', truncG, '.\n')
    
    ## Storage object: matrix with nrow = number of MCMC iterations, and ncol = (1 + p)*truncG, where
    ## truncG the truncation level of the random measure G (an integer given by the values of conc parameter)
    ## (p+1) denoted the number of parameters, each of length truncG, to be stored: p is the number of cluster components (length(tildeVars)) and 1 is for the weights.
    samples <<- matrix(0, nrow = niter, ncol = truncG*(1 + sum(dimTilde))) 
    
    for(iiter in 1:niter){
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
            for(m in 1:dimTilde[j]) {
              uniqueValues[index, jcol] <- values(model, tildeVars[j])[(range[i]-1)*dimTilde[j] + m]
              jcol <- jcol + 1
            }
          }
          index <- index+1
        }
      }
      probs[index] <- concSamples[iiter] 
      newValueIndex <- index 
      
      ## computing G(.) = sum_{l=1}^{truncG} w_l delta_{atom_l} (.):
      # sampling the weights
      vsb <- rbeta(1, 1, concSamples[iiter] + N)
      vsbProd <- 1
      samples[iiter, 1] <<- vsb 
      for(l in 2:(truncG - 1)) {
        vsbProd <- vsbProd * (1-vsb)
        vsb <- rbeta(1, 1, concSamples[iiter]+N)
        samples[iiter, l] <<- vsb * vsbProd 
      }
      samples[iiter, truncG] <<- 1 - sum(samples[iiter, 1:(truncG-1)])
      
      # sampling atoms
      for(l in 1:truncG) {
        ## copy tilde parents into model for use in simulation below when simulate from G_0  
        nimCopy(mvSaved, model, parentNodesTildeVars, row = iiter)
        
        ## first sampled values: w_1 and atom_1
        index <- rcat(prob = probs[1:newValueIndex])
        if(index == newValueIndex){   # sample from G_0
          model$simulate(parentNodesTildeVarsDeps)
          sumDim <- 0
          for(j in 1:p){ 
            for(m in 1:dimTilde[j]) {
              samples[iiter,  truncG + sumDim  + (l-1)*dimTilde[j] + m] <<- values(model, tildeVars[j])[m]
            }
            sumDim <- sumDim + truncG*dimTilde[j]
          }
        } else {   # sample one of the existing values
          jcol <- 1
          sumDim <- 0
          for(j in 1:p){
            for(m in 1:dimTilde[j]) {
              samples[iiter,  truncG + sumDim  + (l-1)*dimTilde[j] + m] <<- uniqueValues[index, jcol] 
              jcol <- jcol + 1
            }
            sumDim <- sumDim + truncG*dimTilde[j]
          }
        }
      }
    }
  },
  methods = list( reset = function () {} )
)



sampleDPmeasure_old <- nimbleFunction(
  name = 'sampleDPmeasure_old',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved){
    
    ## Check if the mvSaved is compiled or not.
    mvIsCompiled <- exists('dll', envir = mvSaved)
    if( mvIsCompiled ) {
      stop("sampleDPmeasure: modelValues object has to be an uncompiled object.\n")
    }
    
    ## Determine variables in the mv object and nodes/variables in the model.
    mvSavedVars <- mvSaved$varNames
    
    stochNodes <- model$getNodeNames(stochOnly = TRUE)
    dataNodes <- model$getNodeNames(dataOnly = TRUE) 
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
    tildeVars <- NULL
    itildeVar <- 1
    dep <- model$getDependencies(targetElements[1], self = FALSE)
    for(i in seq_along(dep)) { 
      expr <- cc_getNodesInExpr(model$getValueExpr(dep[i])) 
      for(j in seq_along(expr)) {
        ## look for cases like thetatilde[xi[i]] to identify 'xi' and extract 'thetaTilde'
        tmpexpr <- parse(text = expr[j])[[1]]
        if(length(tmpexpr) >= 3 && is.call(tmpexpr) && tmpexpr[[1]] == '[') {   
          foundTarget <- all.vars(tmpexpr[[3]]) == dcrpVar   
          if( length(foundTarget) > 0 && sum(foundTarget) > 0 ) {
            tildeVars[itildeVar] <- deparse(tmpexpr[[2]])
            itildeVar <- itildeVar+1 
          }
        }
      }
    }
    if( is.null(tildeVars) )  ## probably unnecessary as checked in CRP sampler, but best to be safe
      stop('sampleDPmeasure: The model should have at least one cluster variable.\n')
    
    ## Check that cluster variables are monitored.
    counts <- tildeVars %in% mvSavedVars
    if( sum(counts) != length(tildeVars) ) 
      stop('sampleDPmeasure: The node(s) representing the cluster variables must be monitored in the MCMC (and therefore stored in the modelValues object).\n')  
    
    ## Check that tilde nodes are continuous and univariate variables (for avoiding long trials in simulating atoms for G in run code:
    if(any(model$isDiscrete(tildeVars)))
      stop('sampleDPmeasure: cluster variables should be continuous random variables.\n')
    
    if(any(model$isMultivariate(tildeVars)))
      stop( 'sampleDPmeasure: only univariate cluster variables are allowed.\n' )
    
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
    
    fixedConc <- TRUE # assume that conc parameter is fixed. This will change in the if statement if necessary
    
    ## Get all parents of xi (membership) variable (i.e., nodes involved in concentration parameter), potentially
    ## needed to determine concentration parameter for each iteration of MCMC
    parentNodesXi <- NULL
    candidateParentNodes <- model$getNodeNames(includeData = FALSE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes == dcrpNode]
    for(i in seq_along(candidateParentNodes)){
      aux <- model$getDependencies(candidateParentNodes[i], self = FALSE) 
      if(sum(aux == dcrpNode)) 
        parentNodesXi <- c(parentNodesXi, candidateParentNodes[i])
    }
    
    if(length(parentNodesXi)) { # concentration parameter is random (or at least a function of other quantities)
      ## Check that values stored in mvSaved are sufficient to calculate dcrp concentration.
      ## Do this by creating a model containing all NAs and then copying in variables that are saved 
      ## and then checking getParam gives a non-NA.
      fixedConc <- FALSE
      
      ## Which parent nodes are saved.
      parentNodesXi <- parentNodesXi[parentNodesXi %in% mvSavedVars]  
      
      ## Create model with NA values
      verbosity <- nimbleOptions('verbose')
      nimbleOptions(verbose = FALSE)
      modelWithNAs <- model$modelDef$newModel(check = FALSE, calculate = FALSE)
      nimbleOptions(verbose = verbosity)
      
      ## Note that this check could fail if current state of model has NAs that result in NA in conc parameter,
      ## but we don't want to use mvSaved as MCMC may not have been run at this point.
      ## Need to think more about this.
      nimCopy(from = model, to = modelWithNAs, nodes = parentNodesXi) 
      if(length(parentNodesXi)) {
        ## These nodes need to be calculated to determine conc param in run code
        parentNodesXiDeps <- model$getDependencies(parentNodesXi, self = FALSE, determOnly = TRUE)  ## excludes dcrpNode
        modelWithNAs$calculate(parentNodesXiDeps)
        dcrpParam <- modelWithNAs$getParam(dcrpNode, 'conc')
        if(is.na(dcrpParam)) 
          stop('sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n') 
      } else stop( 'sampleDPmeasure: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n') 
    } else { ## placeholder since parentNodesXi must only have nodes in the mvSaved for correct compilation
      parentNodesXiDeps <- dcrpNode
      parentNodesXi <- dcrpNode # also used in run code
    }  
    
    N <- length(dataNodes)
    p <- length(tildeVars)
    nTilde <- length(values(model, tildeVars)) / p 
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
    nimCat('sampleDPmeasure: Approximating the random measure by a finite stick-breaking representation with an error smaller than 1e-10, leads to a truncation level of ', truncG, '.\n')
    
    ## Storage object: matrix with nrow = number of MCMC iterations, and ncol = (1 + p)*truncG, where
    ## truncG the truncation level of the random measure G (an integer given by the values of conc parameter)
    ## (p+1) denoted the number of parameters, each of length truncG, to be stored: p is the number of cluster components (length(tildeVars)) and 1 is for the weights.
    samples <<- matrix(0, nrow = niter, ncol = truncG*(p+1)) 
    
    for(iiter in 1:niter){
      ## getting the sampled unique values (tilde variables) and their probabilities of being sampled,
      ## need for computing density later.
      probs <- nimNumeric(N)
      uniqueValues <- matrix(0, nrow = N, ncol = p)  
      xiiter <- mvSaved[dcrpVar, iiter]
      range <- min(xiiter):max(xiiter) 
      index <- 1
      for(i in seq_along(range)){   
        cond <- sum(xiiter == range[i])
        if(cond > 0){
          probs[index] <- cond
          ## slight workaround because can't compile mvSaved[tildeVars[j], iiter]  
          nimCopy(mvSaved, model, tildeVars, row = iiter)
          for(j in 1:p){
            uniqueValues[index, j] <- values(model, tildeVars[j])[range[i]]
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
      paramAux <- numeric(p)
      
      ## copy tilde parents into model for use in simulation below when simulate from G_0  
      nimCopy(mvSaved, model, parentNodesTildeVars, row = iiter)
      
      ## first sampled values: w_1 and atom_1
      index <- rcat(prob = probs[1:newValueIndex])
      if(index == newValueIndex){   # sample from G_0
        model$simulate(parentNodesTildeVarsDeps)
        for(j in 1:p){ 
          samples[iiter, j*truncG + Taux] <<- values(model, tildeVars)[(j-1)*nTilde + 1]
        }
      } else {   # sample one of the existing values
        for(j in 1:p){
          samples[iiter, j*truncG + Taux] <<- uniqueValues[index, j] 
        }
      }
      samples[iiter, Taux] <<- vaux 
      Taux <- Taux + 1
      
      # the rest of the values: w_l and atom_l, l=1, ..truncG-1
      while(Taux <= truncG-1){
        index <- rcat(prob = probs[1:newValueIndex])
        if(index == newValueIndex){  # sample from G_0
          model$simulate(parentNodesTildeVarsDeps)
          for(j in 1:p){ 
            paramAux[j] <- values(model, tildeVars)[(j-1)*nTilde + 1]
          }
        } else{  # sample one of the existing values
          for(j in 1:p){
            paramAux[j] <- uniqueValues[index, j] 
          }
        }
        condaux <- samples[iiter, truncG + 1:(Taux-1)] == paramAux[1]  # check if we sample a new atom or an atom that is in G already
        if(sum(condaux) > 0) { # the atom already exists and we have to update the weights and not include a new value of the params 
          repindex = 1
          while(!condaux[repindex]){
            repindex = repindex + 1
          }
          v1prod <- v1prod * (1-vaux)
          vaux <- rbeta(1, 1, concSamples[iiter]+N)
          samples[iiter, repindex] <<- samples[iiter, repindex] + vaux * v1prod
        } else { # augment the truncation and keep the same parameters
          for(j in 1:p){
            samples[iiter, j*truncG + Taux] <<- paramAux[j]
          }
          v1prod <- v1prod * (1-vaux)
          vaux <- rbeta(1, 1, concSamples[iiter]+N)
          samples[iiter, Taux] <<- vaux * v1prod 
          Taux <- Taux + 1
        }
      }
      
      ## complete the vector of probabilities and atoms: w_Trunc and atom_Trunc
      samples[iiter, truncG] <<- 1 - sum(samples[iiter, 1:(truncG-1)])
      model$simulate(parentNodesTildeVarsDeps)
      for(j in 1:p){ 
        samples[iiter, (j+1)*truncG] <<- values(model, tildeVars)[(j-1)*nTilde + 1]
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
    sample = function(i = integer(), j = integer()) {} ## nothing needed for non-conjugate
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
      model[[marginalizedVar]][j] <<- rnorm(1, postMean, sqrt(postVar)) 
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
      model[[marginalizedVar]][j] <<- rgamma(1, shape = priorShape + y, rate = priorRate + 1)
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
      model[[marginalizedVar]][j] <<- rgamma(1, shape = priorShape + 0.5, rate = priorRate + (y-dataMean)^2/2)
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
      model[[marginalizedVar]][j] <<- rbeta(1, shape1=priorShape1+y, shape2=priorShape2+1-y)
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
      model[[marginalizedVar]][j] <<- rbeta(1, shape1=priorShape1+y, shape2=priorShape2+dataSize-y)
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
      model[[marginalizedVar]][j] <<- rbeta(1, shape1=priorShape1+dataSize, shape2=priorShape2+y)
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
      model[[marginalizedVar]][j] <<- rgamma(1, shape=priorShape+1, rate=priorRate+y)
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
      model[[marginalizedVar]][j] <<- rgamma(1, shape=datashape+priorShape, rate=priorRate+y)
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
      model[[marginalizedVar]][j] <<- rgamma(1, shape=1+priorShape, rate=priorRate+y^dataShape)
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
      model[[marginalizedVar]][j] <<- rgamma(1, shape=dataShape+priorShape, rate=priorRate+1/y)
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
    
    helperFunctions <- nimbleFunctionList(CRP_helper)
    
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
    helperFunctions[[1]] <- eval(as.name(sampler))(model, tildeVars[1], model$expandNodeNames(tildeVars[1]), dataNodes)
    
    curLogProb <- numeric(n)
  },
  
  
  run = function() {
    
    isNonParam <- TRUE # indicates if the MCMC sampling is nonparametric or not. Is nonparametric when there are more cluster parameters than required clusters.
    conc <- model$getParam(target, 'conc')
    helperFunctions[[1]]$storeParams()
    
    xi <- model[[target]]
    
    ## Find unique values in model[[target]]. I'm not relabeling the unique values.
    ## k denotes the number of unique labels
    
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
    while(xiCounts[kNew] > 0 & kNew < n) { # need to make sure don't go beyond length of vector
      kNew <- kNew + 1
    }
    
    
    for(i in 1:n) { # updates one cluster membership at the time , i=1,...,n
      
      xi <- model[[target]]
      xiCounts[xi[i]] <- xiCounts[xi[i]] - 1 # updates counts based on xi[-i]
      
      # computing probability of sampling the k unique values
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
      #  probability of sampling the new label: kNew
      curLogProb[k+1] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i) # probability of sampling a new label
      
      # sampling the index that represents the updated label
      index <- rcat( n=1, exp(curLogProb[1:(k+1)]-max(curLogProb[1:(k+1)])) )
      
      # updating model[[target]][i], xiUniques, xiCounts, kNew
      if(index == (k+1)){ # a new membership variable is sampled
        if( xiCounts[xi[i]] != 0 ) { # a new cluster is created
          if(kNew > min_nTilde) {
            nimCat('CRP_sampler: This MCMC is not fully nonparametric. The MCMC attempted to use more components than the number of cluster parameters.\n')
            kNew <- xi[i]
            isNonParam <- FALSE
          }
          model[[target]][i] <<- kNew 
          k <- k + 1 # a new cluster is created
          xiUniques[k] <- kNew # label of new cluster
          kNew <- kNew + 1
          mySum <- sum(xi == kNew) 
          while(mySum > 0 & kNew < n) { # need to make sure don't go beyond length of vector
            kNew <- kNew+1
            mySum <- sum(xi == kNew)
          }
        } # otherwise, a cluster was deleted so the 'new' label is the deleted label; k, xiUniques, kNew, and model[[target]][i] do not change, update xiCounts
        xiCounts[model[[target]][i]] <- 1
        
        if(isNonParam) { # updating the cluster parameters of the new cluster
          helperFunctions[[1]]$sample(i, model[[target]][i])
        }
      } else { # an existing membership variable is sampled
        model[[target]][i] <<- xiUniques[index]
        xiCounts[xiUniques[index]] <- xiCounts[xiUniques[index]] + 1
        
        if(xiCounts[xi[i]] == 0 ) {  # a cluster is deleted. If no cluster is deleted the two lines above are enough
          ## reorder unique labels because one was deleted
          positionXiInUniques <- which( xiUniques == xi[i] )[1] 
          xiUniques[positionXiInUniques:k] <- xiUniques[(positionXiInUniques+1):(k+1)]
          k <- k - 1
          if(xi[i] < kNew)  ## ensure new cluster is first empty label
            kNew <- xi[i]
        }
      }
    }
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



