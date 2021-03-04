source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)

nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

nimbleProgressBarSetting <- nimbleOptions('MCMCprogressBar')
nimbleOptions(MCMCprogressBar = FALSE)

context('Testing of BNP functionality')



# adding old function of getSamplesDPmeasure function
#--------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------
## old version
getSamplesDPmeasure_old <- function(MCMC, epsilon = 1e-4) {
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
  rsampler <- sampleDPmeasure_old(model, mvSamples, epsilon)
  if(compiled) {
    csampler <- compileNimble(rsampler, project = model)
    csampler$run()
    samplesMeasure <- csampler$samples
  } else {
    rsampler$run()
    samplesMeasure <- rsampler$samples
  }
  
  dcrpVar <- rsampler$dcrpVar
  clusterVarInfo <- nimble:::findClusterNodes(model, dcrpVar) 
  namesVars <- rsampler$tildeVars
  p <- length(namesVars)
  
  truncG <- ncol(samplesMeasure) / (rsampler$tildeVarsColsSum[p+1]+1) 
  namesW <- sapply(seq_len(truncG), function(i) paste0("weight[", i, "]"))
  namesAtoms <- nimble:::getSamplesDPmeasureNames(clusterVarInfo, model, truncG, p)
  
  colnames(samplesMeasure) <- c(namesW, namesAtoms)
  
  output <- list(samples = samplesMeasure, trunc = truncG)
  return(output)
}


sampleDPmeasure_old <- nimbleFunction(
  name = 'sampleDPmeasure_old',

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
    clusterVarInfo <- nimble:::findClusterNodes(model, dcrpVar) 
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
    
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvGammaConjugacy(model, clusterVarInfo, length(dcrpElements), 'dinvgamma'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvGammaConjugacy(model, clusterVarInfo, length(dcrpElements), 'dgamma'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvWishartConjugacy(model, clusterVarInfo, length(dcrpElements), 'dinvwish'))
      isIID <- TRUE
    if(!isIID && length(tildeVars) == 2 && nimble:::checkNormalInvWishartConjugacy(model, clusterVarInfo, length(dcrpElements), 'dwish'))
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
    
    
    ## Storage object to be sized in run code based on MCMC output.
    samples <- matrix(0, nrow = 1, ncol = 1)   
    ## Truncation level of the random measure 
    truncG <- 0 
    niter <- 0
    ## control list extraction
    ## The error of approximation G is given by (conc / (conc +1))^{truncG-1}. 
    ## we are going to define an error of approximation and based on the posterior values of the conc parameter define the truncation level of G
    ## the error is between errors that are considered very very small in the folowing papers
    ## Ishwaran, H., & James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. Journal of the American Statistical Association, 96(453), 161-173.
    ## Ishwaran, H., & Zarepour, M. (2000). Markov chain Monte Carlo in approximate Dirichlet and beta two-parameter process hierarchical models. Biometrika, 87(2), 371-390.
    # epsilon <- 1e-4
    
    setupOutputs(lengthData, dcrpVar)
  },
  
  run=function(){
    
    niter <<- getsize(mvSaved) # number of iterations in the MCMC
    
    # defining the truncation level of the random measure's representation:
    if( fixedConc ) {
      concSamples <- nimNumeric(length = niter, value = model$getParam(dcrpNode, 'conc'))
    } else {
      concSamples <- numeric(niter)
      for( iiter in 1:niter ) {
        nimCopy(from = mvSaved, to = model, nodes = parentNodesXi, row=iiter) 
        model$calculate(parentNodesXiDeps)
        concSamples[iiter] <- model$getParam(dcrpNode, 'conc')
      }
    }
    dcrpAux <- mean(concSamples) + N
    
    truncG <<- log(epsilon) / log(dcrpAux / (dcrpAux+1))
    truncG <<- ceiling(truncG)
    
    ## Storage object: matrix with nrow = number of MCMC iterations, and ncol = (1 + p)*truncG, where
    ## truncG the truncation level of the random measure G (an integer given by the values of conc parameter)
    ## (p+1) denoted the number of parameters, each of length truncG, to be stored: p is the number of cluster components (length(tildeVars)) and 1 is for the weights.
    samples <<- matrix(0, nrow = niter, ncol = truncG*(tildeVarsColsSum[p+1]+1))
    
    ## computing G(.) = sum_{l=1}^{truncG} w_l delta_{atom_l} (.):
    for(iiter in 1:niter){
      checkInterrupt()
      
      ## sampling weights:
      vaux <- rbeta(1, 1, concSamples[iiter] + N)
      v1prod <- 1
      samples[iiter, 1] <<- vaux  
      for(l1 in 2:truncG) {
        v1prod <- v1prod * (1-vaux)
        vaux <- rbeta(1, 1, concSamples[iiter] + N)
        samples[iiter, l1] <<- vaux * v1prod 
      }
      samples[iiter, 1:truncG] <<- samples[iiter, 1:truncG] / (1 - v1prod * (1-vaux)) # normalizing
      
      ## sampling atoms:
      ## getting the sampled unique values (tilde variables) and their probabilities of being sampled,
      ## need for computing density later.
      probs <- nimNumeric(N)
      uniqueValues <- matrix(0, nrow = N, ncol = tildeVarsColsSum[p+1])  
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
            jcols <- (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]
            uniqueValues[index, jcols] <- values(model, tildeVars[j])[mvIndexes[range[i], jcols]]
          }
          index <- index+1
        }
      }
      probs[index] <- concSamples[iiter] 
      newValueIndex <- index 
      
      ## copy tilde parents into model for use in simulation below when simulate atoms of G_0  
      nimCopy(mvSaved, model, parentNodesTildeVars, row = iiter)
      
      sumCol <- truncG
      for(l1 in 1:truncG) {
        index <- rcat(prob = probs[1:newValueIndex])
        if(index == newValueIndex){   # sample from G_0
          model$simulate(parentNodesTildeVarsDeps)
          for(j in 1:p){
            jcols <- (sumCol + 1):(sumCol + tildeVarsCols[j]) 
            samples[iiter, jcols] <<- values(model, tildeVars[j])[mvIndexes[1,  (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]]] # <<-
            sumCol <- sumCol + tildeVarsCols[j] 
          }
        } else {   # sample one of the existing values
          for(j in 1:p){
            jcols <- (sumCol +1):(sumCol + tildeVarsCols[j])
            samples[iiter, jcols] <<- uniqueValues[index, (tildeVarsColsSum[j]+1):tildeVarsColsSum[j+1]] # <<-
            sumCol <- sumCol + tildeVarsCols[j] 
          }
        }
      }
      
    }
  },
  methods = list( reset = function () {} )
)


test_that("Test computations (prior predictive and posterior) and sampler assignment for conjugate CRP samplers", {
  set.seed(0)
  
  # here we test:
  # correct computation of prior predictive and sampling from posterior 
  # correct sampler assignment
  
  ## dnorm_dnorm
  code <- nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        y[i,j] ~ dnorm( mu[xi[i], j] , var = j/2) 
        mu[i, j] ~ dnorm(0.2*j, var=j)
      }
    }
    xi[1:5] ~ dCRP(1, size=5)
  })
  inits <- list(xi = 1:5, 
                mu = matrix(rnorm(5*2, 0), nrow=5,  ncol=2))
  y <- matrix(rnorm(5*2, 10, 1), ncol=2, nrow=5)
  y[4:5, ] <- rnorm(2*2, -10, 1)
  data <- list(y=y)
  model <- nimbleModel(code, data=data, inits=inits,  dimensions=list(mu=c(5,2)), calculate=TRUE)
  mConf <- configureMCMC(model, monitors = c('xi','mu'))  
  mcmc <- buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(c(model$getLogProb('y[1, 1]'), model$getLogProb('y[1, 2]')))
  pT <- sum(c(model$getLogProb('mu[1, 1]'), model$getLogProb('mu[1, 2]')))
  
  dataVar <- c(model$getParam('y[1,1]', 'var') , model$getParam('y[1,2]', 'var') )
  priorVar <- c(model$getParam('mu[1, 1]', 'var'), model$getParam('mu[1, 2]', 'var'))
  priorMean <- c(model$getParam('mu[1, 1]', 'mean') , model$getParam('mu[1, 2]', 'mean'))
  postVar <- 1 / (1 / dataVar + 1 / priorVar) # from conjugate sampler
  postMean <- postVar * (c(data$y[1, 1], data$y[1, 2]) / dataVar + priorMean / priorVar) # from conjugate sampler
  pTgivenY <- dnorm(model$mu[1, 1] , postMean[1], sqrt(postVar[1]), log = TRUE)  + dnorm(model$mu[1, 2] , postMean[2], sqrt(postVar[2]), log = TRUE)# from conjugate sampler
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rnorm(2 , postMean, sqrt(postVar))
  expect_identical(smp, c(model$mu[1, 1], model$mu[1, 2]))
  
  ## conjugate dmnorm_dmnorm 
  code=nimbleCode(
    {
      for(i in 1:5){
        for(j in 1:2) {
          mu[i, 1:2, j] ~ dmnorm(mu0[1:2, j], cov=Cov0[1:2, 1:2, j])
          y[i, 1:2, j] ~ dmnorm(mu[xi[i], 1:2, j], cov=Sigma0[1:2, 1:2, j])  
        }
      }
      xi[1:5] ~ dCRP(conc=1, size=5)
    }
  )
  mu <- array(0, c(5, 2, 2))
  for(j in 1:2) {
    mu[ , ,j] <- matrix(rnorm(5*2, 0, sqrt(0.01)), nrow=5, ncol=2)
  }
  y <- array(0, c(5, 2, 5))
  for(i in 1:5) {
    for(j in 1:2) {
      y[i, ,j] <- rnorm(2, 5, sqrt(0.01))
    }
  }
  mu0 <- matrix(rnorm(2*2), ncol=2, nrow=2)
  Cov0 <- array(0, c(2, 2, 2))
  Sigma0 <- array(0, c(2, 2, 2))
  for(j in 1:2) {
    Cov0[, , j] <- rinvwish_chol(1, chol(matrix(c(10, .7, .7, 10), 2)), 2)
    Sigma0[, , j] <- rinvwish_chol(1, chol(matrix(c(1, .5, .5, 1), 2)), 2)
  }
  
  model = nimbleModel(code, 
                      data = list(y = y),
                      inits = list(xi = 1:5, mu=mu), 
                      constants=list(mu0 =mu0, Cov0 = Cov0, Sigma0 = Sigma0))
  conf <- configureMCMC(model, monitors=c('xi', 'mu'))
  mcmc <- buildMCMC(conf)
  
  # sampler assignment:
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dmnorm_dmnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(c(model$getLogProb('y[1, 1:2, 1]'), model$getLogProb('y[1, 1:2, 2]')))
  pT <- sum(c(model$getLogProb('mu[1, 1:2, 1]'), model$getLogProb('mu[1, 1:2, 2]')))
  
  dataCov <- list(model$getParam('y[1, 1:2, 1]', 'cov') , model$getParam('y[1, 1:2, 2]', 'cov') )
  priorCov <- list(model$getParam('mu[1, 1:2, 1]', 'cov'), model$getParam('mu[1, 1:2, 2]', 'cov'))
  priorMean <- list(model$getParam('mu[1, 1:2, 1]', 'mean') , model$getParam('mu[1, 1:2, 2]', 'mean'))
  
  
  dataPrec <- list(inverse(dataCov[[1]]), inverse(dataCov[[2]]))
  priorPrec <- list(inverse(priorCov[[1]]), inverse(priorCov[[2]]))
  postPrecChol <- list(chol(dataPrec[[1]] + priorPrec[[1]]), chol(dataPrec[[2]] + priorPrec[[2]]))
  postMean <- list(backsolve(postPrecChol[[1]], forwardsolve(t(postPrecChol[[1]]), 
                                                             (dataPrec[[1]] %*% y[1, 1:2, 1] + priorPrec[[1]] %*% priorMean[[1]])[,1])), 
                   backsolve(postPrecChol[[2]], forwardsolve(t(postPrecChol[[2]]), 
                                                             (dataPrec[[2]] %*% y[1, 1:2, 2] + priorPrec[[2]] %*% priorMean[[2]])[,1])))
  pTgivenY <- dmnorm_chol(model$mu[1, 1:2, 1], postMean[[1]], postPrecChol[[1]], prec_param = TRUE, log = TRUE) +
    dmnorm_chol(model$mu[1, 1:2, 2], postMean[[2]], postPrecChol[[2]], prec_param = TRUE, log = TRUE) 
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- rmnorm_chol(1, postMean[[1]], postPrecChol[[1]], prec_param = TRUE) 
  smp2 <- rmnorm_chol(1, postMean[[2]], postPrecChol[[2]], prec_param = TRUE) 
  expect_identical(smp1, model$mu[1, 1:2, 1])
  expect_identical(smp2, model$mu[1, 1:2, 2])
  
  
  ## conjugate dinvgamma_dnorm, 
  code <- nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        y[i,j] ~ dnorm( mu[i, j] , var = s2[xi[i], j]) 
        s2[i, j] ~ dinvgamma(shape = 2*j, scale = 0.1*j) 
      }
    }
    xi[1:5] ~ dCRP(1, size=5)
  })
  inits <- list(xi = 1:5, 
                mu = matrix(rnorm(5*2, 0), nrow=5,  ncol=2),
                s2 = matrix(rinvgamma(5*2, 2, 0.1), nrow=5,  ncol=2))
  y <- matrix(rnorm(5*2, 10, 1), ncol=2, nrow=5)
  y[4:5, ] <- rnorm(2*2, -10, 1)
  data <- list(y=y)
  model <- nimbleModel(code, data=data, inits=inits,  dimensions=list(mu=c(5,2)), calculate=TRUE)
  mConf <- configureMCMC(model, monitors = c('xi','s2'))  
  mcmc <- buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dinvgamma_dnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(c(model$getLogProb('y[1, 1]'), model$getLogProb('y[1, 2]')))
  pT <- sum(c(model$getLogProb('s2[1, 1]'), model$getLogProb('s2[1, 2]')))
  
  dataMean <- c(model$getParam('y[1,1]', 'mean') , model$getParam('y[1,2]', 'mean') )
  priorShape <- c(model$getParam('s2[1, 1]', 'shape'), model$getParam('s2[1, 2]', 'shape'))
  priorScale <- c(model$getParam('s2[1, 1]', 'scale') , model$getParam('s2[1, 2]', 'scale'))
  postShape <- priorShape + 0.5
  postScale <- priorScale + 0.5 * (c(data$y[1, 1], data$y[1, 2]) - dataMean)^2 
  pTgivenY <- dinvgamma(model$s2[1, 1] , shape = postShape[1], scale = postScale[1], log = TRUE)  + 
    dinvgamma(model$s2[1, 2] , shape = postShape[2], scale = postScale[2], log = TRUE) 
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rinvgamma(2 , shape = postShape, scale = postScale)
  expect_identical(smp, c(model$s2[1, 1], model$s2[1, 2]))
  
  ## conjugate dinvwish_dmnorm
  code <- nimbleCode({
    xi[1:5] ~ dCRP(conc = 1, size = 5)
    for(i in 1:5){
      for(j in 1:2) {
        Sigma[1:2, 1:2, i, j] ~ dinvwish(S = R0[1:2, 1:2, j], df = v0[j])
        y[i, 1:2, j] ~ dmnorm(mu[i, 1:2, j],  cov = Sigma[1:2, 1:2, xi[i], j] ) 
      }
    }
  })
  R0 <- array(0, c(2, 2, 2))
  for(j in 1:2) {
    R0[, , j] <- rinvwish_chol(1, chol(matrix(c(10, .7, .7, 10), 2)), 2)
  }
  Sigma <- array(0, c(2,2,5, 2))
  for(i in 1:5){
    for(j in 1:2) {
      Sigma[, , i, j] <- rinvwish_chol(1, chol(matrix(c(1, .5, .5, 1), 2)), 2)
    }
  }
  mu <- array(0, c(5, 2, 2))
  for(j in 1:2) {
    mu[ , ,j] <- matrix(rnorm(5*2, 0, sqrt(0.01)), nrow=5, ncol=2)
  }
  y <- array(0, c(5, 2, 2))
  for(i in 1:5) {
    for(j in 1:2) {
      y[i, ,j] <- rnorm(2, 0, sqrt(0.01))
    }
  }
  data = list(y = y)
  inits = list(xi = 1:5, mu = mu, Sigma = Sigma)
  Consts <- list(v0 = rpois(2, 5), R0 =  R0)
  model = nimbleModel(code, data=data, inits=inits, constants = Consts)
  mConf = configureMCMC(model, monitors = c('xi', 'Sigma'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dinvwish_dmnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  dataMean <- list(model$getParam('y[1, 1:2, 1]', 'mean'), model$getParam('y[1, 1:2, 2]', 'mean'))
  pYgivenT <- sum(model$getLogProb('y[1, 1:2, 1]'), model$getLogProb('y[1, 1:2, 2]'))
  pT <- sum(model$getLogProb('Sigma[1:2, 1:2, 1, 1]'), model$getLogProb('Sigma[1:2, 1:2, 1,  2]'))
  
  df0 <- c(model$getParam('Sigma[1:2, 1:2, 1, 1]', 'df'), model$getParam('Sigma[1:2, 1:2, 1, 2]', 'df'))
  priorScale <- list(model$getParam('Sigma[1:2, 1:2, 1, 1]',  'S'), model$getParam('Sigma[1:2, 1:2, 1, 2]',  'S'))
  
  pTgivenY <- dinvwish_chol(model$Sigma[1:2, 1:2, 1, 1],
                            chol(priorScale[[1]] + (data$y[1, 1:2, 1]-dataMean[[1]])%*%t(data$y[1, 1:2, 1]-dataMean[[1]])),
                            df = (df0[1]+1), scale_param=TRUE, log = TRUE) +
    dinvwish_chol(model$Sigma[1:2, 1:2, 1, 2],
                  chol(priorScale[[2]] + (data$y[1, 1:2, 2]-dataMean[[2]])%*%t(data$y[1, 1:2, 2]-dataMean[[2]])),
                  df = (df0[2]+1), scale_param=TRUE, log = TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)[1]
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- list()
  smp2 <- list()
  smp1[[1]] <- rinvwish_chol(1, chol(priorScale[[1]] + (data$y[1, 1:2, 1]-dataMean[[1]])%*%t(data$y[1, 1:2, 1]-dataMean[[1]])),
                             df = (df0[1]+1), scale_param=TRUE )
  smp1[[2]] <- rinvwish_chol(1, chol(priorScale[[2]] + (data$y[1, 1:2, 2]-dataMean[[2]])%*%t(data$y[1, 1:2, 2]-dataMean[[2]])),
                             df = (df0[2]+1), scale_param=TRUE )
  expect_identical(smp1[[1]], model$Sigma[1:2, 1:2, 1, 1])
  expect_identical(smp1[[2]], model$Sigma[1:2, 1:2, 1, 2])
  
  ## conjugate dwish_dmnorm
  code <- nimbleCode({
    xi[1:5] ~ dCRP(conc = 1, size = 5)
    for(i in 1:5){
      for(j in 1:2) {
        Sigma[1:2, 1:2, i, j] ~ dwish(R = R0[1:2, 1:2, j], df = v0[j])
        y[i, 1:2, j] ~ dmnorm(mu[i, 1:2, j],  prec = Sigma[1:2, 1:2, xi[i], j] ) 
      }
    }
  })
  R0 <- array(0, c(2, 2, 2))
  for(j in 1:2) {
    R0[, , j] <- rwish_chol(1, chol(matrix(c(10, .7, .7, 10), 2)), 2)
  }
  Sigma <- array(0, c(2,2,5, 2))
  for(i in 1:5){
    for(j in 1:2) {
      Sigma[, , i, j] <- rwish_chol(1, chol(matrix(c(1, .5, .5, 1), 2)), 2)
    }
  }
  mu <- array(0, c(5, 2, 2))
  for(j in 1:2) {
    mu[ , ,j] <- matrix(rnorm(5*2, 0, sqrt(0.01)), nrow=5, ncol=2)
  }
  y <- array(0, c(5, 2, 2))
  for(i in 1:5) {
    for(j in 1:2) {
      y[i, ,j] <- rnorm(2, 0, sqrt(0.01))
    }
  }
  data = list(y = y)
  inits = list(xi = 1:5, mu = mu, Sigma = Sigma)
  Consts <- list(v0 = rpois(2, 5), R0 =  R0)
  model = nimbleModel(code, data=data, inits=inits, constants = Consts)
  mConf = configureMCMC(model, monitors = c('xi', 'Sigma'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dwish_dmnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  dataMean <- list(model$getParam('y[1, 1:2, 1]', 'mean'), model$getParam('y[1, 1:2, 2]', 'mean'))
  pYgivenT <- sum(model$getLogProb('y[1, 1:2, 1]'), model$getLogProb('y[1, 1:2, 2]'))
  pT <- sum(model$getLogProb('Sigma[1:2, 1:2, 1, 1]'), model$getLogProb('Sigma[1:2, 1:2, 1,  2]'))
  
  df0 <- c(model$getParam('Sigma[1:2, 1:2, 1, 1]', 'df'), model$getParam('Sigma[1:2, 1:2, 1, 2]', 'df'))
  priorScale <- list(model$getParam('Sigma[1:2, 1:2, 1, 1]',  'R'), model$getParam('Sigma[1:2, 1:2, 1, 2]',  'R'))
  
  pTgivenY <- dwish_chol(model$Sigma[1:2, 1:2, 1, 1],
                         chol(priorScale[[1]] + (data$y[1, 1:2, 1]-dataMean[[1]])%*%t(data$y[1, 1:2, 1]-dataMean[[1]])),
                         df = (df0[1]+1), scale_param=FALSE, log = TRUE) +
    dwish_chol(model$Sigma[1:2, 1:2, 1, 2],
               chol(priorScale[[2]] + (data$y[1, 1:2, 2]-dataMean[[2]])%*%t(data$y[1, 1:2, 2]-dataMean[[2]])),
               df = (df0[2]+1), scale_param=FALSE, log = TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)[1]
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- list()
  smp2 <- list()
  smp1[[1]] <- rwish_chol(1, chol(priorScale[[1]] + (data$y[1, 1:2, 1]-dataMean[[1]])%*%t(data$y[1, 1:2, 1]-dataMean[[1]])),
                          df = (df0[1]+1), scale_param=FALSE )
  smp1[[2]] <- rwish_chol(1, chol(priorScale[[2]] + (data$y[1, 1:2, 2]-dataMean[[2]])%*%t(data$y[1, 1:2, 2]-dataMean[[2]])),
                          df = (df0[2]+1), scale_param=FALSE )
  expect_identical(smp1[[1]], model$Sigma[1:2, 1:2, 1, 1])
  expect_identical(smp1[[2]], model$Sigma[1:2, 1:2, 1, 2])
  
  ## conjugate dnorm_invgamma_dnorm 
  code = nimbleCode({
    xi[1:5] ~ dCRP(conc=1, size=5)
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dnorm(j, var = s2[i, j]/kappa[j])
        s2[i, j] ~ dinvgamma(shape=j+1, scale=j)
        y[i, j] ~ dnorm(mu[xi[i], j], var=s2[xi[i], j])  
      }
    }
    for(j in 1:2) {
      kappa[j] <- 2+j
    }
  })
  y <- matrix(rnorm(5*2, 10, 1), ncol=2, nrow=5)
  y[4:5, ] <- rnorm(2*2, -10, 1)
  data = list(y = y)
  inits = list(xi = 1:5, mu=matrix(rnorm(5*2), ncol=2, nrow=5), s2=matrix(rinvgamma(5*2, 2, 1), ncol=2, nrow=5))
  model = nimbleModel(code, data=data, inits=inits)
  mConf = configureMCMC(model, monitors = c('xi','mu', 's2'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(model$getLogProb('y[1, 1]'), model$getLogProb('y[1, 2]'))
  pT1 <- sum(model$getLogProb('mu[1,1]'),  model$getLogProb('mu[1, 2]'))
  pT2 <- sum(model$getLogProb('s2[1, 1]'), model$getLogProb('s2[1, 2]'))
  
  priorMean <- c(model$getParam('mu[1, 1]', 'mean'), model$getParam('mu[1, 2]', 'mean')) 
  kappa <- c(values(model, 's2[1, 1]')[1]/model$getParam('mu[1, 1]', 'var'), values(model, 's2[1, 2]')[1]/model$getParam('mu[1, 2]', 'var'))
  priorShape <- c(model$getParam('s2[1, 1]', 'shape'), model$getParam('s2[1, 2]', 'shape'))
  priorScale <- c(model$getParam('s2[1, 1]',  'scale'), model$getParam('s2[1, 2]',  'scale'))
  pTgivenY2 <- dinvgamma(model$s2[1, 1], shape = priorShape[1] + 1/2,
                         scale = priorScale[1] + kappa[1] * (data$y[1,1] - priorMean[1])^2 / (2*(1+kappa[1])),
                         log=TRUE) + dinvgamma(model$s2[1, 2], shape = priorShape[2] + 1/2,
                                               scale = priorScale[2] + kappa[2] * (data$y[1,2] - priorMean[2])^2 / (2*(1+kappa[2])),
                                               log=TRUE)
  pTgivenY1 <- dnorm(model$mu[1, 1], mean = (kappa[1] * priorMean[1] + data$y[1, 1])/(1 + kappa[1]), 
                     sd = sqrt(model$s2[1, 1] / (1+kappa[1])),
                     log=TRUE)  + dnorm(model$mu[1, 2], mean = (kappa[2] * priorMean[2] + data$y[1, 2])/(1 + kappa[2]), 
                                        sd = sqrt(model$s2[1, 2] / (1+kappa[2])),
                                        log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT1 + pT2 + pYgivenT - pTgivenY1 - pTgivenY2)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- rep(0,2) # to preserve order in sampling
  smp2 <- rep(0,2)
  smp1[1] <- rinvgamma(1, shape = priorShape[1] + 1/2,
                       scale = priorScale[1] + kappa[1] * (data$y[1,1] - priorMean[1])^2 / (2*(1+kappa[1])))
  smp2[1] <- rnorm(1, mean = (kappa[1] * priorMean[1] + data$y[1,1])/(1 + kappa[1]), 
                   sd = sqrt(smp1[1] / (1+kappa[1])))
  smp1[2] <- rinvgamma(1, shape = priorShape[2] + 1/2,
                       scale = priorScale[2] + kappa[2] * (data$y[1,2] - priorMean[2])^2 / (2*(1+kappa[2])) )
  smp2[2] <- rnorm(1, mean = (kappa[2] * priorMean[2] + data$y[1, 2])/(1 + kappa[2]), 
                   sd = sqrt(smp1[2] / (1+kappa[2]))) 
  expect_identical(smp1, c(model$s2[1, 1], model$s2[1, 2]))
  expect_identical(smp2, c(model$mu[1, 1], model$mu[1, 2]) )
  
  ## conjugate dnorm_gamma_dnorm model
  code = nimbleCode({
    xi[1:5] ~ dCRP(conc=1, size=5)
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dnorm(j, tau = s2[i, j]*kappa[j])
        s2[i, j] ~ dgamma(shape=j+1, rate=j)
        y[i, j] ~ dnorm(mu[xi[i], j], tau=s2[xi[i], j])  
      }
    }
    for(j in 1:2) {
      kappa[j] <- 2+j
    }
  })
  y <- matrix(rnorm(5*2, 10, 1), ncol=2, nrow=5)
  y[4:5, ] <- rnorm(2*2, -10, 1)
  data = list(y = y)
  inits = list(xi = 1:5, mu=matrix(rnorm(5*2), ncol=2, nrow=5), s2=matrix(rgamma(5*2, 1, 2), ncol=2, nrow=5))
  model = nimbleModel(code, data=data, inits=inits)
  mConf = configureMCMC(model, monitors = c('xi','mu', 's2'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_gamma_dnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(model$getLogProb('y[1, 1]'), model$getLogProb('y[1, 2]'))
  pT1 <- sum(model$getLogProb('mu[1,1]'),  model$getLogProb('mu[1, 2]'))
  pT2 <- sum(model$getLogProb('s2[1, 1]'), model$getLogProb('s2[1, 2]'))
  
  priorMean <- c(model$getParam('mu[1, 1]', 'mean'), model$getParam('mu[1, 2]', 'mean')) 
  kappa <- c(model$getParam('mu[1, 1]', 'tau') / values(model, 's2[1, 1]')[1], model$getParam('mu[1, 2]', 'tau') / values(model, 's2[1, 2]')[1])
  priorShape <- c(model$getParam('s2[1, 1]', 'shape'), model$getParam('s2[1, 2]', 'shape'))
  priorRate <- c(model$getParam('s2[1, 1]',  'rate'), model$getParam('s2[1, 2]',  'rate'))
  pTgivenY2 <- dgamma(model$s2[1, 1], shape = priorShape[1] + 1/2,
                      rate = priorRate[1] + kappa[1] * (data$y[1,1] - priorMean[1])^2 / (2*(1+kappa[1])),
                      log=TRUE) + dgamma(model$s2[1, 2], shape = priorShape[2] + 1/2,
                                         rate = priorRate[2] + kappa[2] * (data$y[1,2] - priorMean[2])^2 / (2*(1+kappa[2])),
                                         log=TRUE)
  pTgivenY1 <- dnorm(model$mu[1, 1], mean = (kappa[1] * priorMean[1] + data$y[1, 1])/(1 + kappa[1]), 
                     sd = sqrt(1/(model$s2[1, 1] *(1+kappa[1]))),
                     log=TRUE)  + dnorm(model$mu[1, 2], mean = (kappa[2] * priorMean[2] + data$y[1, 2])/(1 + kappa[2]), 
                                        sd = sqrt(1/(model$s2[1, 2] *(1+kappa[2]))),
                                        log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT1 + pT2 + pYgivenT - pTgivenY1 - pTgivenY2)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- rep(0,2) # to preserve order in sampling
  smp2 <- rep(0,2)
  smp1[1] <- rgamma(1, shape = priorShape[1] + 1/2,
                    rate = priorRate[1] + kappa[1] * (data$y[1,1] - priorMean[1])^2 / (2*(1+kappa[1])))
  smp2[1] <- rnorm(1, mean = (kappa[1] * priorMean[1] + data$y[1,1])/(1 + kappa[1]), 
                   sd = sqrt(1 / (smp1[1]*(1+kappa[1]))))
  smp1[2] <- rgamma(1, shape = priorShape[2] + 1/2,
                    rate = priorRate[2] + kappa[2] * (data$y[1,2] - priorMean[2])^2 / (2*(1+kappa[2])) )
  smp2[2] <- rnorm(1, mean = (kappa[2] * priorMean[2] + data$y[1, 2])/(1 + kappa[2]), 
                   sd = sqrt(1 / (smp1[2]*(1+kappa[2])))) 
  expect_identical(smp1, c(model$s2[1, 1], model$s2[1, 2]))
  expect_identical(smp2, c(model$mu[1, 1], model$mu[1, 2]) )
  
  ## conjugate dmnorm_invwish_dmnorm
  code <- nimbleCode({
    xi[1:5] ~ dCRP(conc = 1, size = 5)
    for(i in 1:5){
      for(j in 1:2) {
        Sigma[1:2, 1:2, i, j] ~ dinvwish(S = R0[1:2, 1:2, j], df = v0[j])
        SigmaAux[1:2, 1:2, i, j] <- Sigma[1:2, 1:2, i, j]  / k0[j]
        
        mu[i, 1:2, j] ~ dmnorm(mu0[1:2, j], cov = SigmaAux[1:2, 1:2, i, j] )
        y[i, 1:2, j] ~ dmnorm(mu[xi[i], 1:2, j],  cov = Sigma[1:2, 1:2, xi[i], j] )  
      }
    }
  })
  R0 <- array(0, c(2, 2, 2))
  for(j in 1:2) {
    R0[, , j] <- rinvwish_chol(1, chol(matrix(c(10, .7, .7, 10), 2)), 2)
  }
  Sigma <- array(0, c(2,2,5, 2))
  for(i in 1:5){
    for(j in 1:2) {
      Sigma[, , i, j] <- rinvwish_chol(1, chol(matrix(c(1, .5, .5, 1), 2)), 2)
    }
  }
  mu <- array(0, c(5, 2, 2))
  for(j in 1:2) {
    mu[ , ,j] <- matrix(rnorm(5*2, 0, sqrt(0.01)), nrow=5, ncol=2)
  }
  y <- array(0, c(5, 2, 2))
  for(i in 1:5) {
    for(j in 1:2) {
      y[i, ,j] <- rnorm(2, 0, sqrt(0.01))
    }
  }
  data = list(y = y)
  inits = list(xi = 1:5, mu = mu, Sigma = Sigma)
  Consts <- list(mu0 = matrix(rnorm(4), ncol=2, nrow=2), v0 = rpois(2, 5),
                 k0 = rgamma(2, 1, 1),
                 R0 =  R0)
  model = nimbleModel(code, data=data, inits=inits, constants = Consts)
  mConf = configureMCMC(model, monitors = c('xi','mu', 'Sigma'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dmnorm_invwish_dmnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(model$getLogProb('y[1, 1:2, 1]'), model$getLogProb('y[1, 1:2, 2]'))
  pT1 <- sum(model$getLogProb('mu[1,1:2, 1]'),  model$getLogProb('mu[1, 1:2, 2]'))
  pT2 <- sum(model$getLogProb('Sigma[1:2, 1:2, 1, 1]'), model$getLogProb('Sigma[1:2, 1:2, 1,  2]'))
  
  priorMean <- list(model$getParam('mu[1, 1:2, 1]', 'mean'), model$getParam('mu[1, 1:2, 2]', 'mean')) 
  kappa <- c(values(model, 'Sigma[1:2, 1:2, 1, 1]')[1]/model$getParam('mu[1, 1:2, 1]', 'cov')[1, 1], 
             values(model, 'Sigma[1:2, 1:2, 1, 2]')[1]/model$getParam('mu[1, 1:2, 2]', 'cov')[1, 1])
  df0 <- c(model$getParam('Sigma[1:2, 1:2, 1, 1]', 'df'), model$getParam('Sigma[1:2, 1:2, 1, 2]', 'df'))
  priorScale <- list(model$getParam('Sigma[1:2, 1:2, 1, 1]',  'S'), model$getParam('Sigma[1:2, 1:2, 1, 2]',  'S'))
  
  pTgivenY2 <- dinvwish_chol(model$Sigma[1:2, 1:2, 1, 1],
                             chol(priorScale[[1]] + (kappa[1]/(kappa[1]+1)) * (data$y[1, 1:2, 1]-priorMean[[1]])%*%t(data$y[1, 1:2, 1]-priorMean[[1]])),
                             df = (df0[1]+1), scale_param=TRUE, log = TRUE) +
    dinvwish_chol(model$Sigma[1:2, 1:2, 1, 2],
                  chol(priorScale[[2]] + (kappa[2]/(kappa[2]+1)) * (data$y[1, 1:2, 2]-priorMean[[2]])%*%t(data$y[1, 1:2, 2]-priorMean[[2]])),
                  df = (df0[2]+1), scale_param=TRUE, log = TRUE)
  pTgivenY1 <- dmnorm_chol(model$mu[1, 1:2, 1], mean = (kappa[1] * priorMean[[1]] + data$y[1, 1:2, 1])/(1 + kappa[1]), 
                           chol( model$Sigma[1:2, 1:2, 1, 1] / (1+kappa[1]) ),
                           prec_param = FALSE, log = TRUE) + 
    dmnorm_chol(model$mu[1, 1:2, 2], mean = (kappa[2] * priorMean[[2]] + data$y[1, 1:2, 2])/(1 + kappa[2]), 
                chol( model$Sigma[1:2, 1:2, 1, 2] / (1+kappa[2]) ),
                prec_param = FALSE, log = TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)[1]
  
  expect_equal(pY, pT1 + pT2 + pYgivenT - pTgivenY1 - pTgivenY2)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- list()
  smp2 <- list()
  smp1[[1]] <- rinvwish_chol(1, chol(priorScale[[1]] + (kappa[1]/(kappa[1]+1)) * (data$y[1, 1:2, 1]-priorMean[[1]])%*%t(data$y[1, 1:2, 1]-priorMean[[1]])),
                             df = (df0[1]+1), scale_param=TRUE )
  smp2[[1]] <- rmnorm_chol(1, mean = (kappa[1] * priorMean[[1]] + data$y[1, 1:2, 1])/(1 + kappa[1]), 
                           chol( smp1[[1]] / (1+kappa[1]) ), prec_param = FALSE)
  smp1[[2]] <- rinvwish_chol(1, chol(priorScale[[2]] + (kappa[2]/(kappa[2]+1)) * (data$y[1, 1:2, 2]-priorMean[[2]])%*%t(data$y[1, 1:2, 2]-priorMean[[2]])),
                             df = (df0[2]+1), scale_param=TRUE )
  smp2[[2]] <- rmnorm_chol(1, mean = (kappa[2] * priorMean[[2]] + data$y[1, 1:2, 2])/(1 + kappa[2]), 
                           chol( smp1[[2]] / (1+kappa[2]) ), prec_param = FALSE)
  expect_identical(smp1[[1]], model$Sigma[1:2, 1:2, 1, 1])
  expect_identical(smp1[[2]], model$Sigma[1:2, 1:2, 1, 2])
  expect_identical(smp2[[1]], model$mu[1, 1:2, 1])
  expect_identical(smp2[[2]], model$mu[1, 1:2, 2])
  
  ## conjugate dmnorm_wish_dmnorm
  code <- nimbleCode({
    xi[1:5] ~ dCRP(conc = 1, size = 5)
    for(i in 1:5){
      for(j in 1:2) {
        Sigma[1:2, 1:2, i, j] ~ dwish(R = R0[1:2, 1:2, j], df = v0[j])
        SigmaAux[1:2, 1:2, i, j] <- Sigma[1:2, 1:2, i, j]  * k0[j]
        
        mu[i, 1:2, j] ~ dmnorm(mu0[1:2, j], prec = SigmaAux[1:2, 1:2, i, j] )
        y[i, 1:2, j] ~ dmnorm(mu[xi[i], 1:2, j],  prec = Sigma[1:2, 1:2, xi[i], j] )  
      }
    }
  })
  R0 <- array(0, c(2, 2, 2))
  for(j in 1:2) {
    R0[, , j] <- rwish_chol(1, chol(matrix(c(10, .7, .7, 10), 2)), 2)
  }
  Sigma <- array(0, c(2,2,5, 2))
  for(i in 1:5){
    for(j in 1:2) {
      Sigma[, , i, j] <- rwish_chol(1, chol(matrix(c(1, .5, .5, 1), 2)), 2)
    }
  }
  mu <- array(0, c(5, 2, 2))
  for(j in 1:2) {
    mu[ , ,j] <- matrix(rnorm(5*2, 0, sqrt(0.01)), nrow=5, ncol=2)
  }
  y <- array(0, c(5, 2, 2))
  for(i in 1:5) {
    for(j in 1:2) {
      y[i, ,j] <- rnorm(2, 0, sqrt(0.01))
    }
  }
  data = list(y = y)
  inits = list(xi = 1:5, mu = mu, Sigma = Sigma)
  Consts <- list(mu0 = matrix(rnorm(4), ncol=2, nrow=2), v0 = rpois(2, 5),
                 k0 = rgamma(2, 1, 1),
                 R0 =  R0)
  model = nimbleModel(code, data=data, inits=inits, constants = Consts)
  mConf = configureMCMC(model, monitors = c('xi','mu', 'Sigma'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dmnorm_wish_dmnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(model$getLogProb('y[1, 1:2, 1]'), model$getLogProb('y[1, 1:2, 2]'))
  pT1 <- sum(model$getLogProb('mu[1,1:2, 1]'),  model$getLogProb('mu[1, 1:2, 2]'))
  pT2 <- sum(model$getLogProb('Sigma[1:2, 1:2, 1, 1]'), model$getLogProb('Sigma[1:2, 1:2, 1,  2]'))
  
  priorMean <- list(model$getParam('mu[1, 1:2, 1]', 'mean'), model$getParam('mu[1, 1:2, 2]', 'mean')) 
  kappa <- c(model$getParam('mu[1, 1:2, 1]', 'prec')[1, 1]/values(model, 'Sigma[1:2, 1:2, 1, 1]')[1], 
             model$getParam('mu[1, 1:2, 2]', 'prec')[1, 1]/values(model, 'Sigma[1:2, 1:2, 1, 2]')[1])
  df0 <- c(model$getParam('Sigma[1:2, 1:2, 1, 1]', 'df'), model$getParam('Sigma[1:2, 1:2, 1, 2]', 'df'))
  priorRate <- list(model$getParam('Sigma[1:2, 1:2, 1, 1]',  'R'), model$getParam('Sigma[1:2, 1:2, 1, 2]',  'R'))
  
  pTgivenY2 <- dwish_chol(model$Sigma[1:2, 1:2, 1, 1],
                          chol(priorRate[[1]] + (kappa[1]/(kappa[1]+1)) * (data$y[1, 1:2, 1]-priorMean[[1]])%*%t(data$y[1, 1:2, 1]-priorMean[[1]])),
                          df = (df0[1]+1), scale_param=FALSE, log = TRUE) +
    dwish_chol(model$Sigma[1:2, 1:2, 1, 2],
               chol(priorRate[[2]] + (kappa[2]/(kappa[2]+1)) * (data$y[1, 1:2, 2]-priorMean[[2]])%*%t(data$y[1, 1:2, 2]-priorMean[[2]])),
               df = (df0[2]+1), scale_param=FALSE, log = TRUE)
  pTgivenY1 <- dmnorm_chol(model$mu[1, 1:2, 1], mean = (kappa[1] * priorMean[[1]] + data$y[1, 1:2, 1])/(1 + kappa[1]), 
                           chol( model$Sigma[1:2, 1:2, 1, 1] * (1+kappa[1]) ),
                           prec_param = TRUE, log = TRUE) + 
    dmnorm_chol(model$mu[1, 1:2, 2], mean = (kappa[2] * priorMean[[2]] + data$y[1, 1:2, 2])/(1 + kappa[2]), 
                chol( model$Sigma[1:2, 1:2, 1, 2] * (1+kappa[2]) ),
                prec_param = TRUE, log = TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)[1]
  
  expect_equal(pY, pT1 + pT2 + pYgivenT - pTgivenY1 - pTgivenY2)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp1 <- list()
  smp2 <- list()
  smp1[[1]] <- rwish_chol(1, chol(priorRate[[1]] + (kappa[1]/(kappa[1]+1)) * (data$y[1, 1:2, 1]-priorMean[[1]])%*%t(data$y[1, 1:2, 1]-priorMean[[1]])),
                          df = (df0[1]+1), scale_param=FALSE )
  smp2[[1]] <- rmnorm_chol(1, mean = (kappa[1] * priorMean[[1]] + data$y[1, 1:2, 1])/(1 + kappa[1]), 
                           chol( smp1[[1]] * (1+kappa[1]) ), prec_param = TRUE)
  smp1[[2]] <- rwish_chol(1, chol(priorRate[[2]] + (kappa[2]/(kappa[2]+1)) * (data$y[1, 1:2, 2]-priorMean[[2]])%*%t(data$y[1, 1:2, 2]-priorMean[[2]])),
                          df = (df0[2]+1), scale_param=FALSE )
  smp2[[2]] <- rmnorm_chol(1, mean = (kappa[2] * priorMean[[2]] + data$y[1, 1:2, 2])/(1 + kappa[2]), 
                           chol( smp1[[2]] * (1+kappa[2]) ), prec_param = TRUE)
  expect_identical(smp1[[1]], model$Sigma[1:2, 1:2, 1, 1])
  expect_identical(smp1[[2]], model$Sigma[1:2, 1:2, 1, 2])
  expect_identical(smp2[[1]], model$mu[1, 1:2, 1])
  expect_identical(smp2[[2]], model$mu[1, 1:2, 2])
  
  ## conjugate dgamma_dnorm, 
  code <- nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        y[i,j] ~ dnorm( mu[i, j] , tau = s2[xi[i], j]) 
        s2[i, j] ~ dgamma(shape = j, rate = j+1) 
      }
    }
    xi[1:5] ~ dCRP(1, size=5)
  })
  inits <- list(xi = 1:5, 
                mu = matrix(rnorm(5*2, 0), nrow=5,  ncol=2),
                s2 = matrix(rgamma(5*2, 0.1, rate=1), nrow=5,  ncol=2))
  y <- matrix(rnorm(5*2, 10, 1), ncol=2, nrow=5)
  data <- list(y=y)
  model <- nimbleModel(code, data=data, inits=inits,  dimensions=list(mu=c(5,2)), calculate=TRUE)
  mConf <- configureMCMC(model, monitors = c('xi','s2'))  
  mcmc <- buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dnorm")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(c(model$getLogProb('y[1, 1]'), model$getLogProb('y[1, 2]')))
  pT <- sum(c(model$getLogProb('s2[1, 1]'), model$getLogProb('s2[1, 2]')))
  
  dataMean <- c(model$getParam('y[1,1]', 'mean') , model$getParam('y[1,2]', 'mean') )
  priorShape <- c(model$getParam('s2[1, 1]', 'shape'), model$getParam('s2[1, 2]', 'shape'))
  priorRate <- c(model$getParam('s2[1, 1]', 'rate') , model$getParam('s2[1, 2]', 'rate'))
  postShape <- priorShape + 0.5
  postRate <- priorRate + 0.5 * (c(data$y[1, 1], data$y[1, 2]) - dataMean)^2 
  pTgivenY <- dgamma(model$s2[1, 1] , shape = postShape[1], rate = postRate[1], log = TRUE)  + 
    dgamma(model$s2[1, 2] , shape = postShape[2], rate = postRate[2], log = TRUE) 
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- rgamma(2 , shape = postShape, rate = postRate)
  expect_identical(smp, c(model$s2[1, 1], model$s2[1, 2]))
  
  ## dbeta_dbern
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dbeta(1+j,j)
        y[i, j] ~ dbern(mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y = matrix(rbinom(10, size=1, prob=0.1), ncol=2, nrow=5)
  y[4:5, ] <- rbinom(4, size=1, prob=0.9)
  data = list(y=y)
  inits = list(xi = 1:5, mu=matrix(rbeta(10, 1, 1), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits = inits)
  mConf = configureMCMC(m, monitors = c('xi','mu'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dbern")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1]'), m$getLogProb('y[1, 2]'))
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape1 <- c(m$getParam('mu[1, 1]', 'shape1'), m$getParam('mu[1, 2]', 'shape1'))
  priorShape2 <- c(m$getParam('mu[1, 1]', 'shape2'), m$getParam('mu[1, 2]', 'shape2'))
  pTgivenY <- dbeta(m$mu[1, 1], shape1=priorShape1[1]+data$y[1, 1], shape2=priorShape2[1]+1-data$y[1, 1], log=TRUE) +
    dbeta(m$mu[1, 2], shape1=priorShape1[2]+data$y[1, 2], shape2=priorShape2[2]+1-data$y[1, 2], log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rbeta(1 , shape1=priorShape1[1]+data$y[1, 1], shape2=priorShape2[1]+1-data$y[1, 1]),
           rbeta(1 , shape1=priorShape1[2]+data$y[1, 2], shape2=priorShape2[2]+1-data$y[1, 2]) )
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]))
  
  ## dbeta_dbin
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dbeta(j,5+j)
        y[i, j] ~ dbinom(size=10, prob=mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y = matrix(rbinom(10, size=10, prob=0.1), ncol=2, nrow=5)
  data = list(y=y)
  inits = list(xi = 1:5, mu=matrix(rbeta(10, 1, 1), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits = inits)
  mConf = configureMCMC(m, monitors = c('xi','mu'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dbin")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1]'), m$getLogProb('y[1, 2]'))
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape1 <- c(m$getParam('mu[1, 1]', 'shape1'), m$getParam('mu[1, 2]', 'shape1'))
  priorShape2 <- c(m$getParam('mu[1, 1]', 'shape2'), m$getParam('mu[1, 2]', 'shape2'))
  dataSize <- c(m$getParam('y[1, 1]', 'size'), m$getParam('y[1, 2]', 'size'))
  pTgivenY <- dbeta(m$mu[1, 1], shape1=priorShape1[1]+data$y[1, 1], shape2=priorShape2[1]+dataSize[1]-data$y[1, 1], log=TRUE)+
    dbeta(m$mu[1, 2], shape1=priorShape1[2]+data$y[1, 2], shape2=priorShape2[2]+dataSize[2]-data$y[1, 2], log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rbeta(1 , shape1=priorShape1[1]+data$y[1, 1], shape2=priorShape2[1]+dataSize[1]-data$y[1, 1]),
           rbeta(1 , shape1=priorShape1[2]+data$y[1, 2], shape2=priorShape2[2]+dataSize[2]-data$y[1, 2]) )
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]))
  
  ## dbeta_dnegbin
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dbeta(j,j+1)
        y[i, j] ~ dnegbin(size=10, prob=mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y = matrix(rnbinom(10, size=10, prob=0.1), ncol=2, nrow=5)
  data = list(y=y)
  inits = list(xi = 1:5, mu=matrix(rbeta(10, 1, 1), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits= inits)
  mConf = configureMCMC(m, monitors = c('xi','mu'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dbeta_dnegbin")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1]'), m$getLogProb('y[1, 2]'))
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape1 <- c(m$getParam('mu[1, 1]', 'shape1'), m$getParam('mu[1, 2]', 'shape1'))
  priorShape2 <- c(m$getParam('mu[1, 1]', 'shape2'), m$getParam('mu[1, 2]', 'shape2'))
  dataSize <- c(m$getParam('y[1, 1]', 'size'), m$getParam('y[1, 2]', 'size'))
  pTgivenY <- dbeta(m$mu[1, 1], shape1=priorShape1[1]+dataSize[1], shape2=priorShape2[1]+data$y[1, 1], log=TRUE) +
    dbeta(m$mu[1, 2], shape1=priorShape1[2]+dataSize[2], shape2=priorShape2[2]+data$y[1, 2], log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rbeta(1 , shape1=priorShape1[1]+dataSize[1], shape2=priorShape2[1]+data$y[1, 1]), 
           rbeta(1 , shape1=priorShape1[2]+dataSize[2], shape2=priorShape2[2]+data$y[1, 2]))
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]) )
  

  ## dgamma_dpois
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dgamma(j,j+1)
        y[i, j] ~ dpois(mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y = matrix(rpois(10, 1), ncol=2, nrow=5)
  data = list(y=y)
  inits = list(xi = 1:5, mu=matrix(rgamma(10, 1, 5), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits= inits)
  cm<-compileNimble(m)
  mConf = configureMCMC(m, monitors=c('mu', 'xi'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dpois")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1]'), m$getLogProb('y[1, 2]'))
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape <- c(m$getParam('mu[1, 1]', 'shape'), m$getParam('mu[1, 2]', 'shape'))
  priorRate <- c(m$getParam('mu[1, 1]', 'rate'), m$getParam('mu[1, 2]', 'rate'))
  pTgivenY <- dgamma(m$mu[1, 1], shape = priorShape[1] + data$y[1, 1], rate = priorRate[1] + 1, log=TRUE) +
    dgamma(m$mu[1, 2], shape = priorShape[2] + data$y[1, 2], rate = priorRate[2] + 1, log=TRUE) 
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rgamma(1 , shape = priorShape[1] + data$y[1, 1], rate = priorRate[1] + 1), 
           rgamma(1 , shape = priorShape[2] + data$y[1, 2], rate = priorRate[2] + 1))
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]))

 
  ## dgamma_dexp:
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dgamma(j,j+1)
        y[i, j] ~ dexp(mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y = matrix(rexp(10, 1), ncol=2, nrow=5)
  data = list(y=y)
  inits = list(xi = 1:5, mu=matrix(rgamma(10, 1, 1), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits= inits)
  mConf = configureMCMC(m, monitors = c('xi','mu'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dexp")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1]'),m$getLogProb('y[1, 2]')) 
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape <- c(m$getParam('mu[1, 1]', 'shape'), m$getParam('mu[1, 2]', 'shape'))
  priorRate <- c(m$getParam('mu[1, 1]', 'rate'), m$getParam('mu[1, 2]', 'rate'))
  pTgivenY <- dgamma(m$mu[1,1], shape=priorShape[1]+1, rate=priorRate[1]+data$y[1, 1], log=TRUE)+
    dgamma(m$mu[1, 2], shape=priorShape[2]+1, rate=priorRate[2]+data$y[1, 2], log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rgamma(1, shape=priorShape[1]+1, rate=priorRate[1]+data$y[1, 1]), 
           rgamma(1, shape=priorShape[2]+1, rate=priorRate[2]+data$y[1, 2]))
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]))
  
  ## dgamma_dgamma:
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dgamma(j, rate = j+1)
        y[i, j] ~ dgamma(4, rate = mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y = matrix(rgamma(10, 4, 4), ncol=2, nrow=5)
  data = list(y = y)
  inits = list(xi = 1:5, mu=matrix(rgamma(10, 1, 5), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits= inits)
  mConf = configureMCMC(m, monitors = c('xi','mu'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dgamma")
  
  pYgivenT <- sum(m$getLogProb('y[1, 1]'), m$getLogProb('y[1, 2]'))
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape <- c(m$getParam('mu[1, 1]', 'shape'), m$getParam('mu[1, 2]', 'shape'))
  priorRate <- c(m$getParam('mu[1, 1]', 'rate'), m$getParam('mu[1, 2]', 'rate'))
  dataShape <- c(m$getParam('y[1, 1]', 'shape'), m$getParam('y[1, 2]', 'shape'))
  pTgivenY <- dgamma(m$mu[1, 1], shape=dataShape[1]+priorShape[1], rate=priorRate[1]+data$y[1, 1], log=TRUE) +
    dgamma(m$mu[1, 2], shape=dataShape[2]+priorShape[2], rate=priorRate[2]+data$y[1, 2], log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rgamma(1, shape=dataShape[1]+priorShape[1], rate=priorRate[1]+data$y[1, 1]), 
           rgamma(1, shape=dataShape[2]+priorShape[2], rate=priorRate[2]+data$y[1, 2]))
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]))
  
  ## dgamma_dweib:
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dgamma(j, 5+j)
        y[i, j] ~ dweib(shape=4*j, lambda = mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y <- matrix(rweibull(10, 4, 4), ncol=2, nrow=5)
  data = list(y = y)
  inits = list(xi = 1:5, mu=matrix(rgamma(10, 1, 5), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits= inits)
  mConf = configureMCMC(m, monitors=list('xi', 'mu'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dweib")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1]'), m$getLogProb('y[1, 2]'))
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape <- c(m$getParam('mu[1, 1]', 'shape'), m$getParam('mu[1, 2]', 'shape'))
  priorRate <- c(m$getParam('mu[1, 1]', 'rate'), m$getParam('mu[1, 2]', 'rate'))
  dataShape <- c(m$getParam('y[1, 1]', 'shape'), m$getParam('y[1, 2]', 'shape'))
  pTgivenY <- dgamma(m$mu[1, 1], shape=1+priorShape[1], rate=priorRate[1]+data$y[1,1]^dataShape[1], log=TRUE) +
    dgamma(m$mu[1, 2], shape=1+priorShape[2], rate=priorRate[2]+data$y[1, 2]^dataShape[2], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rgamma(1, shape=1+priorShape[1], rate=priorRate[1]+data$y[1, 1]^dataShape[1]),
           rgamma(1, shape=1+priorShape[2], rate=priorRate[2]+data$y[1, 2]^dataShape[2]))
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]))
  
  
  ## dgamma_dinvgamma:
  code = nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        mu[i, j] ~ dgamma(j, rate=5+j)
        y[i, j] ~ dinvgamma(shape=4*j, scale = mu[xi[i], j]) 
      }
    }
    xi[1:5] ~ dCRP(conc=1, size=5)
  })
  y <- matrix(rinvgamma(10, 4, 3), ncol=2, nrow=5)
  data = list(y = y)
  inits = list(xi = 1:5, mu=matrix(rgamma(10, 1, 5), ncol=2, nrow=5))
  m = nimbleModel(code, data=data, inits= inits)
  mConf = configureMCMC(m, monitors = list('xi', 'mu'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dinvgamma")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1]'), m$getLogProb('y[1, 2]'))
  pT <- sum(m$getLogProb('mu[1, 1]'), m$getLogProb('mu[1, 2]'))
  
  priorShape <- c(m$getParam('mu[1, 1]', 'shape'), m$getParam('mu[1, 2]', 'shape'))
  priorRate <- c(m$getParam('mu[1, 1]', 'rate'), m$getParam('mu[1, 2]', 'rate'))
  dataShape <- c(m$getParam('y[1, 1]', 'shape'), m$getParam('y[1, 2]', 'shape'))
  pTgivenY <- dgamma(m$mu[1, 1], shape=dataShape[1]+priorShape[1], rate=priorRate[1]+1/data$y[1, 1], log=TRUE)+
    dgamma(m$mu[1, 2], shape=dataShape[2]+priorShape[2], rate=priorRate[2]+1/data$y[1, 2], log=TRUE)
  
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- c(rgamma(1, shape=dataShape[1]+priorShape[1], rate=priorRate[1]+1/data$y[1, 1]), 
           rgamma(1, shape=dataShape[2]+priorShape[2], rate=priorRate[2]+1/data$y[1, 2]))
  expect_identical(smp, c(m$mu[1, 1], m$mu[1, 2]))
  
  
  ## ddirch_dmulti:
  code=nimbleCode(
    {
      for(i in 1:5){
        for(j in 1:2) {
          p[i, 1:3, j] ~ ddirch(alpha=alpha0[1:3, j])
          y[i, 1:3, j] ~ dmulti(prob=p[xi[i], 1:3, j], size=3)  
        }
      }
      xi[1:5] ~ dCRP(conc=1, size=5)
    }
  )
  alpha0 <- matrix(rgamma(3*2, 1, 1), ncol=2, nrow=3)
  p <- array(0, c(5, 3, 2))
  for(i in 1:5) {
    for(j in 1:2) {
      p[i, , j] <- rdirch(1, c(1, 1, 1))
    }
  }
  
  y <- array(0, c(5, 3, 2))
  for(i in 1:5){
    for(j in 1:2) {
      y[i, , j] = rmulti(1, prob=c(0.01,0.01,0.98), size=3) 
    }
  }
  data = list(y = y)
  m = nimbleModel(code, 
                  data = data,
                  inits = list(xi = 1:5, p=p), 
                  constants=list(alpha0 = alpha0))
  mConf = configureMCMC(m, monitors = list('xi', 'p'))
  mcmc = buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  # computation of prior predictive and posterior sampling of parameters:
  pYgivenT <- sum(m$getLogProb('y[1, 1:3, 1]'), m$getLogProb('y[1, 1:3, 2]'))
  pT <- sum(m$getLogProb('p[1, 1:3, 1]'), m$getLogProb('p[1, 1:3, 2]'))
  
  priorAlpha <- list(m$getParam('p[1, 1:3, 1]', 'alpha'), m$getParam('p[1, 1:3, 2]', 'alpha'))
  pTgivenY <- ddirch(m$p[1,1:3, 1], alpha = priorAlpha[[1]]+data$y[1, 1:3, 1], log=TRUE) + 
    ddirch(m$p[1,1:3, 2], alpha = priorAlpha[[2]]+data$y[1, 1:3, 2], log=TRUE)
  
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$storeParams()
  pY <- mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)
  
  expect_equal(pY, pT + pYgivenT - pTgivenY)
  
  set.seed(1)
  mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]]$sample(1, 1)
  set.seed(1)
  smp <- list(rdirch(1, alpha = priorAlpha[[1]]+data$y[1, 1:3, 1]),
              rdirch(1, alpha = priorAlpha[[2]]+data$y[1, 1:3, 2]))
  expect_identical(smp[[1]], m$p[1, 1:3, 1])
  expect_identical(smp[[2]], m$p[1, 1:3, 2])
  
  
 
  }
)
  
test_that("sampleDPmeasure: testing that required variables in MCMC modelValues are monitored", {
  set.seed(1)
  
  ## membership variable not being monitored
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 ~ dgamma(1, 1)
    for(i in 1:6){
      mu[i] ~ dnorm(0, 1)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6, conc0 = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The node having the dCRP distribution')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'conc0', 'mu'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_silent(output <- getSamplesDPmeasure(mMCMC))
  
  ## cluster variable not being monitored
  code <- nimbleCode({
    xi[1:6] ~ dCRP(1, 6)
    mu0 ~ dnorm(0, 1)
    s20 ~ dgamma(1, 1)
    for(i in 1:6){
      mu[i] ~ dnorm(mu0, s20)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6, mu0 = 0, s20 = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The node\\(s\\) representing the cluster variables') 
  
  mConf <- configureMCMC(m, monitors = c('mu', 'xi'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes')
  
  mConf <- configureMCMC(m, monitors = c('mu', 'xi', 'mu0', 's20'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_silent(output <- getSamplesDPmeasure(mMCMC))
  
  ## concentration parameter not being monitored:
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 ~ dgamma(a, rate=b)
    a ~ dgamma(1, rate=1)
    b ~ dgamma(1, rate=0.1)
    for(i in 1:6){
      mu[i] ~ dnorm(0, 1)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6, conc0 = 1, a = 1, b = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m, monitors = c('xi', 'mu'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'conc0'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  outputG <- getSamplesDPmeasure(mMCMC)
  
  ## concentration parameter deterministic parent not being monitored:
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 <- a + b
    a ~ dgamma(1, rate=1)
    b <- d + 1
    d ~ dgamma(1, 1)
    for(i in 1:6){
      mu[i] ~ dnorm(0, 1)
      y[i] ~ dnorm(mu[xi[i]], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2), mu = 1:6,  a = 1,d=1,  conc0=1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m, monitors = c('xi', 'mu'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'conc0'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'b'))
  mMCMC <- buildMCMC(mConf)
  expect_error(getSamplesDPmeasure(mMCMC),
               'sampleDPmeasure: The stochastic parent nodes of the membership')
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'b', 'd'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_silent(outputG <- getSamplesDPmeasure(mMCMC))
  
  mConf <- configureMCMC(m, monitors = c('xi', 'mu', 'a', 'd'))
  mMCMC <- buildMCMC(mConf)
  expect_message(output <- runMCMC(mMCMC, niter=1))
  expect_silent(outputG <- getSamplesDPmeasure(mMCMC))
})


test_that("check iid assumption in sampleDPmeasure", {
  set.seed(1)
  
  ## univariate cluster parameters are not iid
  code <- nimbleCode({
    for(i in 1:10){
      muTilde[i] ~ dnorm(i, 1)
      y[i] ~ dnorm(muTilde[xi[i]], 1)
    }
    xi[1:10] ~ dCRP(conc = 1, size=10)
  })
  Inits <- list( xi = sample(1:2, size=10, replace=TRUE), 
                 muTilde = rep(1, 10))
  Data <- list(y = c(rnorm(10, 0,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors = c('muTilde','xi'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m, showCompilerOutput = FALSE)
  output <- runMCMC(cMCMC,  niter=1, nburnin = 0, thin=1)
  expect_error(samplesG <- getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }

  ## also not IID
  code=nimbleCode({
      xi[1:3] ~ dCRP(1, size = 3)
      thetatilde[1] ~ dnorm(0, 1)
      thetatilde[2] ~ dt(0, 1, 1)
      thetatilde[3] ~ dt(0, 1, 1)
      s2tilde[1] ~ dinvgamma(2, 1)
      s2tilde[2] ~ dgamma(1, 1)
      s2tilde[3] ~ dgamma(1, 1)
      for(i in 1:3){
          y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])
      }
  }
  )
  Inits <- list(xi = rep(1, 3), thetatilde=rep(0,3), s2tilde=rep(1,3))
  Data <- list(y = rnorm(3,-5, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi'))
  expect_silent(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  output <- cMCMC$run(1)
  expect_error(getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }

  
  ## one cluster param not with same distribution
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      thetatilde[1] ~ dnorm(0, 1)
      thetatilde[2] ~ dt(0, 1, 1)
      thetatilde[3] ~ dt(0, 1, 1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), thetatilde=rep(0,3))
  Data=list(y=rnorm(10, 0,1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  cMCMC$run(1)
  expect_error(getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }

  
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      thetatilde[1] ~ dnorm(0, 1)
      thetatilde[2] ~ dt(0, 1, 1)
      thetatilde[3] ~ dt(0, 1, 1)
      s2tilde[1] ~ dinvgamma(2, 1)
      s2tilde[2] ~ dgamma(1, 1)
      s2tilde[3] ~ dgamma(1, 1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])
      }
    }
  )
  Inits=list(xi=rep(1, 10), thetatilde=rep(0,3), s2tilde=rep(1,3))#
  Data=list(y=rnorm(10, 0,1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde', 's2tilde', 'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  cMCMC$run(1, reset=FALSE) 
  expect_error(getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
  
  # bivariate cluster parameters are not iid case in a model with multiple observations per cluster ID
  code=nimbleCode(
    {
      for(j in 1:3) {
        for(i in 1:4){
          muj[j, 1:2, i] <- (i+j)*mu0[1:2]
          muTilde[j, 1:2, i] ~ dmnorm(muj[j, 1:2, i], cov=Cov0[1:2, 1:2])
          y[j, 1:2, i] ~ dmnorm(muTilde[j, 1:2, xi[i]], cov=Sigma0[1:2, 1:2])  
        }
      }
      xi[1:4] ~ dCRP(conc=1, size=4)
    }
  )
  muTilde <- array(0, c(3, 2, 4))
  for(j in 1:3) {
    muTilde[ j, ,] <- matrix(0, nrow=4, ncol=2)
  }
  y <- array(0, c(3, 2, 4))
  for(i in 1:2) {
    for(j in 1:2) {
      y[j, ,i] <- rnorm(2, 5, sqrt(0.01))
    }
    y[3, ,i] <- rnorm(2,10, sqrt(0.01))
  }
  for(i in 3:4) {
    for(j in 1:2) {
      y[j, ,i] <- rnorm(2, -5, sqrt(0.01))
    }
    y[3, ,i] <- rnorm(2, -10, sqrt(0.01))
  }
  m = nimbleModel(code, 
                  data = list(y = y),
                  inits = list(xi = 1:4, muTilde=muTilde), 
                  constants=list(mu0 = rep(0,2), Cov0 = diag(10, 2), Sigma0 = diag(1, 2)))
  cmodel <- compileNimble(m)
  conf <- configureMCMC(m, monitors=c('xi', 'muTilde'))
  mcmc <- buildMCMC(conf)
  cMCMC <- compileNimble(mcmc, project = m)
  cMCMC$run(1)
  expect_error(getSamplesDPmeasure(cMCMC),
               'sampleDPmeasure: cluster parameters have to be independent and identically')
  if(.Platform$OS.type != "windows") {
    nimble:::clearCompiled(m)
  }
  
})


test_that("check use of epsilon parameter in getSamplesDPmeasure", {
  set.seed(1)
  
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
    sd0 ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)      
    mu0 ~ dnorm(0, var=10)
  })
  
  n <- 30
  constants <- list(n = n)
  data <- list(y = rnorm(n, 0, 1))
  inits <- list(alpha = 1, mu0 = 0, sd0 = 5, xi = 1:n,
                muTilde = rep(0,n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  
  outputG <- getSamplesDPmeasure(cmcmc, setSeed = 1)
  tr1 <- nrow(outputG[[1]])
  
  outputG <- getSamplesDPmeasure(cmcmc, epsilon = 0.1, setSeed = 1)
  tr2 <- nrow(outputG[[1]])
  
  outputG <- getSamplesDPmeasure(cmcmc, epsilon = 0.00001, setSeed = 1)
  tr3 <- nrow(outputG[[1]])

  expect_true(tr1 > tr2,
              info='getSamplesDPmeasure: truncation level for larger epsilon incorrectly computed')
  expect_true(tr1 < tr3,
              info='getSamplesDPmeasure: truncation level for smaller epsilon incorrectly computed')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }
})





test_that("Test opening of new clusters in CRP sampler ", {

  ## test for updating new cluster parameters, with conjugacy
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. 
  inits <- list(alpha = 1, mu0 = 0, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)
  
  ## now check that cmodel$muTilde[2] is near 50 and first obs has moved to 2nd cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_equal(output[1, 'muTilde[2]'], 50, tolerance = 2, info = 'incorrect update of parameter for second cluster',
               check.attributes = FALSE)
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=2), 'incorrect cluster for first obs')
  expect_identical(output[1, 'muTilde[1]'], c('muTilde[1]'=0), 'incorrect update of parameter for first cluster')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }

  ## test for updating new cluster parameters, without conjugacy
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. prior generates better value for first data point so new cluster should be opened.
  inits <- list(alpha = 1, mu0 = 50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] has changed and first obs has moved to 2nd cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_true(output[1, 'muTilde[2]'] != -50, 'incorrect update of parameter for second cluster')
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=2), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }

  ## test for (not) updating new cluster parameters, without conjugacy, no movement expected
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. prior generates even worse value for new cluster in terms of first obs, so no movement expected.
  inits <- list(alpha = 1, mu0 = -50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, 50, rep(0, n-2)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] is unchanged and first obs has stayed in only cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_true(output[1, 'muTilde[2]'] == 50, 'incorrect update of parameter for second cluster')
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=1), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }
  
  ## test for updating new cluster parameters, without conjugacy, movement to existing cluster expected
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. Prior generates even worse value for new cluster in terms of first obs, so no movement expected.
  inits <- list(alpha = 1, mu0 = -50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  inits$xi[1] <- 2
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] has remained and first obs has moved to first cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_true(output[1, 'muTilde[2]'] == -50, 'incorrect update of parameter for second cluster')
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=1), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }

  ## test for updating new cluster parameters, without conjugacy, singleton with movement to new cluster, same label, expected
  code <- nimbleCode({
    for(i in 1:n) {
      y[i] ~ T(dnorm(mu[i], 1), -500, 500) # force non-conjugacy
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
    }
    xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 20
  constants <- list(n = n)
  ## all data plausibly from first cluster except 1st data point
  data <- list(y = c(50, rep(0, n-1)))
  ## muTilde is good for all but first data point. prior generates better value for new cluster in terms of first obs, so movement expected, but singleton, so new cluster should be same ID as old cluster.
  inits <- list(alpha = 1, mu0 = 50, sd0 = 5, xi = rep(1, n),
                muTilde = c(0, -50, rep(0, n-2)))
  inits$xi[1] <- 2
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cmodel <- compileNimble(model)
  conf <- configureMCMC(model, monitors = c('xi', 'muTilde', 'sd0', 'alpha', 'mu0'))
  conf$removeSamplers('muTilde')
  mcmc <- buildMCMC(conf)
  cmcmc <- compileNimble(mcmc, project = model)

  ## now check that cmodel$muTilde[2] has changed and first obs has 'stayed' in 2nd cluster
  set.seed(1)
  output <- runMCMC(cmcmc, niter=1, nburnin=0, thin=1 , inits=inits, setSeed=FALSE)
  expect_equal(output[1, 'muTilde[2]'], 50, tolerance = 5,
               info = 'incorrect update of parameter for second cluster',
               check.attributes = FALSE)
  expect_identical(output[1, 'xi[1]'], c('xi[1]'=2), 'incorrect cluster for first obs')
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(model)
  }
})


test_that("Test reset frunction in CRP sampler ", {
  set.seed(1)
  
  ## the data set has more clusters than cluster parameters are available
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:2)
        thetatilde[i] ~ dnorm(0, 1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), thetatilde=rep(0,2))#
  Data=list(y=c(rnorm(3,-5, 1), rnorm(3,5, 1), rnorm(4, 0,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  expect_output(cMCMC$run(1), info='CRP_sampler: This MCMC is for a parametric model')
  cMCMC$run(1, reset=FALSE)
  if(.Platform$OS.type != "windows") {
      nimble:::clearCompiled(m)
  }
})    

test_that("Test that not nonparametric MCMC message in CRP sampler is printed", {
  set.seed(1)
  
  ## 2 mean cluster parameters, data is mixture of 3 normals, concentration param is fixed
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:2)
        thetatilde[i] ~ dnorm(mean=0, var=10)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), 
             thetatilde=c(0,0))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m)
  expect_output(out <- runMCMC(mcmc=cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is for a parametric model.')
  
  
  ## 2 mean cluster parameters, data is mixture of 3 normals, concentration param is random
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(conc0 , size=10)
      conc0 ~ dgamma(1, 1)
      for(i in 1:2)
        thetatilde[i] ~ dnorm(mean=0, var=10)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1, 10), 
             thetatilde=c(0,0), conc0=1)
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m)
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m)
  expect_output(out <- runMCMC(mcmc=cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is not for a proper model.')
  

  ## fewer cluster parameters than observations, concentration param is fixed, conjugate case
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:5)
        thetatilde[i] ~ dnorm(mean=0, var=10)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1:5, 2), 
             thetatilde=rep(0,5))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  expect_output(out <- runMCMC(cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is for a parametric model.')
  
  
  ## fewer cluster parameters than onservation, concentration param is fixed, non conjugate case
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:5)
        thetatilde[i] ~ dt(0,1,1)
      for(i in 1:10){
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=rep(1:5, 2), 
             thetatilde=rep(0,5))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  expect_warning(mMCMC <- buildMCMC(mConf))
  cMCMC <- compileNimble(mMCMC, project = m) 
  expect_output(out <- runMCMC(cMCMC, niter=1, nburnin = 0, thin=1),
                'CRP_sampler: This MCMC is for a parametric model.')
  
  
   
  ## no message is sent when xi=1:n
  code=nimbleCode(
    {
      xi[1:10] ~ dCRP(1 , size=10)
      for(i in 1:10){
        thetatilde[i] ~ dnorm(mean=0, var=10)
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  Inits=list(xi=1:10, 
             thetatilde=rep(0,10))
  Data=list(y=c(rnorm(3,-5, 1), rnorm(4, 0, 1), rnorm(3,5,1)))
  m <- nimbleModel(code, data=Data, inits=Inits)
  cm <- compileNimble(m)
  mConf <- configureMCMC(m, monitors =  c('thetatilde',  'xi'))
  mMCMC <- buildMCMC(mConf)
  cMCMC <- compileNimble(mMCMC, project = m)
  expect_silent(out <- runMCMC(cMCMC, niter=1, nburnin = 0, thin=1))  
  
  
  ## dirichlet-multinomial model, no message is sent when xi = 1:n
  code=nimbleCode(
    {
      for(i in 1:4){
        p[i,1:3] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[xi[i],1:3], size=3)
      }
      xi[1:4] ~ dCRP(conc=1, size=4)
    }
  )
  p0 <- matrix(0, ncol=3, nrow=4)
  y0 <- matrix(0, ncol=3, nrow=4)
  for(i in 1:4){
    p0[i,]=rdirch(1, c(1, 1, 1))
    y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
  }
  m = nimbleModel(code, 
                  data = list(y = y0),
                  inits = list(xi = 1:4, p=p0), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m, monitors=c('p', 'xi'))
  mcmc <- buildMCMC(conf)
  cm = compileNimble(m)
  cmcmc=compileNimble(mcmc,project=m)
  expect_silent(cmcmc$run(100))
  
})


test_that("Check error given when model has no cluster variables", {
    ## Originally this tested whether there are no tildeVars but with new check for 'xi' appearing
    ## in non-index role, the error is caught differently.
  set.seed(1)
  code <- nimbleCode({
    xi[1:6] ~ dCRP(conc0, 6)
    conc0 ~ dgamma(1, 1)
    for(i in 1:6){
      y[i] ~ dnorm(xi[i], 1)
    }
  })
  Inits <- list(xi = c(1,1,2,1,1,2),  conc0 = 1)
  Data <- list( y =  rnorm(6))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  
  expect_error(buildMCMC(mConf) ,
               'sampler_CRP: Detected that the CRP variable is used in some way not as an index')
  
})


test_that("dCRP nimble function calculates density correctly",{
  
  x <- c(1,1,2,1,1,2)
  conc <- 1
  
  truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
    (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1))
  ltruth <- log(truth)
  
  expect_equal(dCRP(x, conc, size=length(x), log=FALSE),
               truth,
               info = paste0("incorrect dCRP nimble function calculation"))
  
  expect_equal(dCRP(x, conc, size=length(x), log=TRUE),
               ltruth,
               info = paste0("incorrect dCRP nimble function calculation in log scale"))
  
  cdCRP <- compileNimble(dCRP)
  
  expect_equal(cdCRP(x, conc, size=length(x)), (truth), 
               info = paste0("incorrect dCRP value in compiled nimble function"))
  
  expect_equal(cdCRP(x, conc, size=length(x), log=TRUE), (ltruth), 
               info = paste0("incorrect dCRP value in compiled nimble function  in log scale"))
  
  expect_equal(dCRP(x, conc=-1, size=length(x), log=FALSE),
               NaN,
               info = paste0("incorrect parameters space allowed"))
  
  expect_error(dCRP(x, conc=1, size=3, log=FALSE), "length of 'x' has to be equal to 'size'")
  
  expect_error(dCRP(x, conc=1, size=10, log=FALSE), "length of 'x' has to be equal to 'size'")
  
})


test_that("CRP model calculation and dimensions are correct:", {
  
  x <- c(1,1,2,1,1,2)
  conc <- 1
  
  truth <- (conc/(conc+1-1))*(1/(conc+2-1))*(conc/(conc+3-1))*
    (2/(conc+4-1))*(3/(conc+5-1))*(1/(conc+6-1))
  ltruth <- log(truth)
  
  CRP_code <- nimbleCode({
    x[1:6] ~ dCRP(conc, size=6)
  })
  
  Consts <- list(conc = 1)
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model <- nimbleModel(CRP_code, data=Inits, constants=Consts)
  
  CRP_model$x <- x
  expect_equal(exp(CRP_model$calculate()), truth,
               info = paste0("incorrect likelihood value for dCRP"))
  
  c_CRP_model <- compileNimble(CRP_model)
  c_CRP_model$x
  expect_equal(exp(c_CRP_model$calculate()), truth,
               info = paste0("incorrect likelihood value for compiled dCRP"))
  
  
  ## different length of x and size:
  CRP_code2 <- nimbleCode({
    x[1:6] ~ dCRP(1, size=10)
  })
  
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model2 <- nimbleModel(CRP_code2, data=Inits)
  expect_error(CRP_model2$calculate(), "length of 'x' has to be equal to 'size'")
  
  ## different length of x and size:
  CRP_code3 <- nimbleCode({
    x[1:6] ~ dCRP(1, size=3)
  })
  
  Inits <- list(x = c(1,1,2,1,1,2))
  CRP_model3 <- nimbleModel(CRP_code3, data=Inits)
  expect_error(CRP_model3$calculate(), "length of 'x' has to be equal to 'size'")
    
})


test_that("random sampling from CRP in model with additional levels", {
  
  conc <- 1
  set.seed(0)
  size <- 6
  r_samps <- t(replicate(10000, rCRP(n = 1, conc, size = size)))
  ## K is the number of unique components in x of length 6
  true_EK <- sum(conc/(conc+1:size-1))
  
  expect_equal(mean(apply(r_samps, 1, function(x)length(unique(x)))), true_EK, 
               tol = 0.01,
               info = "Difference in expected mean of K exceeds tolerance")
  
  ## sampling from the model:
  set.seed(1)
  CRP_code <- nimbleCode({
    x[1:6] ~ dCRP(conc=1, size=6)
    for(i in 1:6){
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[x[i]], 1)
    }
  })
  Inits <- list(x = c(1,1,2,1,1,2), mu = 1:6)
  Data <- list( y =  rnorm(6))
  CRP_model <- nimbleModel(CRP_code, data=Data, inits=Inits)
  c_CRP_model <- compileNimble(CRP_model)
  
  simul_samp <- function(model) {
    model$simulate()
    return(model$x)
  }
  simul_samps <- t(replicate(10000, simul_samp(c_CRP_model)))
  
  expect_equal(mean(apply(simul_samps, 1, function(x)length(unique(x)))), true_EK, 
               tol = 0.01,
               info = "Difference in expected mean of K, from compiled model, exceeds tolerance")
  
})


# these tests need to be updated: these are based on CRP sampler and now we are using CRP_moreGeneral
test_that("Testing conjugacy detection with models using CRP", { 
  
  ## dnorm_dnorm with truncation
  code = nimbleCode({
    for(i in 1:4) 
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    for(i in 1:2){
      mu[i] ~ dnorm(0,1)}
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "sampler_CRP: The number of clusters based on the cluster parameters is less")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm one more level of hierarchy
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(beta,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
    beta ~ dnorm(0,1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4), beta =1))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  
  ## dnorm_dnorm and deterministic nodes
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dnorm(mui[i], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm and deterministic nodes and truncation
  code = nimbleCode({
    for(i in 1:4) {
      mui[i] <- mu[xi[i]]
      y[i] ~ dnorm(mui[i], sd = 1)
    }
    for(i in 1:2){
      mu[i] ~ dnorm(0,1)}
    xi[1:4] ~ dCRP(conc=1, size=4) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm and non-standard indexing
  code = nimbleCode({
    for(i in 1:4) {
      mu[i, 2] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i], 2], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=cbind(rnorm(4),rnorm(4))))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dnorm and non-standard indexing
  code = nimbleCode({
    for(i in 1:4) {
      mu[2, i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[2, xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=t(cbind(rnorm(4),rnorm(4)))))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  
  ## dnorm_dpois
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dpois(10)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rpois(4, 10)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## dnorm_invgamma; we skip conjugacy if data nodes have scaled variance
  code = nimbleCode({
    for(i in 1:4) {
      s2tilde[i] ~ dinvgamma(a,b)
      s2[i] <- lambda * s2tilde[xi[i]]
      y[i] ~ dnorm(0, var = s2[i])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
    lambda ~ dgamma(1, 1)
    a ~ dgamma(1, 1)
    b ~ dgamma(1, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), s2=rinvgamma(4, 1,1), a=1, b=1, lambda=2))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## conjugate normal-normal model mu(i,j) are not iid: non conjugate sampler is assigned
  code <- nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        y[i,j] ~ dnorm( mu[xi[i], j] , var = 1) 
        mu[i, j] ~ dnorm(i+j, var=100) # thetaTilde_{i,1:J} are iid
      }
    }
    xi[1:5] ~ dCRP(1, size=5)
  })
  inits <- list(xi = rep(1, 5), 
                mu = matrix(rnorm(5*2, 0), nrow=5,  ncol=2))
  y <- matrix(rnorm(5*2, 10, 1), ncol=2, nrow=5)
  y[4:5, ] <- rnorm(2*2, -10, 1)
  data <- list(y=y)
  model <- nimbleModel(code, data=data, inits=inits,  dimensions=list(mu=c(5,2)), calculate=TRUE)
  mConf <- configureMCMC(model, monitors = c('xi','mu'))  
  mcmc <- buildMCMC(mConf)
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## non conjugate normal-invgamma-normal model, mu(i,j) and sigma2(i,j) iid cluster params
  code <- nimbleCode({
    for(i in 1:5) {
      for(j in 1:2) {
        y[i,j] ~ dnorm( mu[xi[i], j] , var = sigma2[xi[i], j]) 
        mu[i, j] ~ dnorm(0, var=100) #  iid
        sigma2[i, j] ~ dinvgamma(2, 1) #  iid
      }
    }
    xi[1:5] ~ dCRP(1, size=5)
  })
  inits <- list(xi = rep(1, 5), 
                mu = matrix(rnorm(5*2), nrow=5,  ncol=2), 
                sigma2 = matrix(rinvgamma(5*2, 2, 1), nrow=5,  ncol=2))
  y <- matrix(rnorm(5*2, 10, 1), ncol=2, nrow=5)
  y[4:5, ] <- rnorm(2*2, -10, 1)
  data <- list(y=y)
  model <- nimbleModel(code, data=data, inits=inits,  dimensions=list(mu=c(5,2), sigma2=c(5,2)), calculate=TRUE)
  cmodel<-compileNimble(model)
  mConf <- configureMCMC(model, monitors = c('xi','mu', 'sigma2'))  
  mcmc <- buildMCMC(mConf)
  
  # sampler assignment:
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
  
  ## non-standard ordering/indexing of ddirch-dmulti
  code=nimbleCode(
    {
      for(i in 1:4){
        p[1:3, i] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[1:3, xi[i]], size=3)
      }
      xi[1:4] ~ dCRP(conc=1, size=4)
    }
  )
  set.seed(1)
  p0 <- matrix(0, ncol=3, nrow=4)
  y0 <- matrix(0, ncol=3, nrow=4)
  for(i in 1:4){
    p0[i,]=rdirch(1, c(1, 1, 1))
    y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
  }
  m = nimbleModel(code, 
                  data = list(y = y0),
                  inits = list(xi = rep(1,4), p=t(p0)), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc <- buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  code=nimbleCode(
    {
      for(i in 1:4){
        p[i, 2:4] ~ ddirch(alpha=alpha0[1:3])
        y[i,1:3] ~ dmulti(prob=p[xi[i], 2:4], size=3)
      }
      xi[1:4] ~ dCRP(conc=1, size=4)
    }
  )
  set.seed(1)
  p0 <- matrix(0, ncol=3, nrow=4)
  y0 <- matrix(0, ncol=3, nrow=4)
  for(i in 1:4){
    p0[i,]=rdirch(1, c(1, 1, 1))
    y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
  }
  p0 <- cbind(rep(0, 4), p0)
  m = nimbleModel(code, 
                  data = list(y = y0),
                  inits = list(xi = rep(1,4), p=p0), 
                  constants=list(alpha0 = c(1,1,1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc <- buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_ddirch_dmulti")
  
  ## dnorm, dinvgamma, not conjugate
  code = nimbleCode({
    for(i in 1:4) {
      s2[i] ~ dinvgamma(1, 1)
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], var = s2[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4), s2=rinvgamma(4, 1,1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## dnorm, dinvgamma, not conjugate
  code = nimbleCode({
    for(i in 1:4) {
      sigma[i] ~ dinvgamma(1, 1)
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = sigma[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), sigma = rinvgamma(4, 1,1)))
  conf <- configureMCMC(m)
  mcmc=buildMCMC(conf)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## dnorm_invgamma, conjugate; we detect conjugacy
  code = nimbleCode({
    for(i in 1:4) {
      s2[i] ~ dinvgamma(a,b)
      mu[i] ~ dnorm(0, var = s2[i]/kappa)
      y[i] ~ dnorm(mu[xi[i]], var = s2[xi[i]])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
    kappa ~ dgamma(1, 1)
    a ~ dgamma(1, 1)
    b ~ dgamma(1, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu=rnorm(4), s2=rinvgamma(4, 1,1), a=1, b=1, kappa=2))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(nimble:::checkCRPconjugacy(m, 'xi[1:4]'), "conjugate_dnorm_invgamma_dnorm")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  
  ## model with deterministic nodes
  code = nimbleCode({
    for(i in 1:4) {
      s2Tilde[i] ~ dinvgamma(a,b)
      s2[i] <- s2Tilde[xi[i]]
      muTilde[i] ~ dnorm(0, var = s2Tilde[i]/kappa)
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], var = s2[i])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
    kappa ~ dgamma(1, 1)
    a ~ dgamma(1, 1)
    b ~ dgamma(1, 1)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), muTilde=rnorm(4), s2Tilde=rinvgamma(4, 1,1), a=1, b=1, kappa=2))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(nimble:::checkCRPconjugacy(m, 'xi[1:4]'), "conjugate_dnorm_invgamma_dnorm")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  
  
  ## dgamma_dexp and deterministic nodes
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dexp(mui[i])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dgamma_dexp")
  
  
  ## dgamma_dexp, deterministic nodes, and conjugacy is broken
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dexp(mui[i]+3)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
  ## dgamma_dexp, deterministic nodes, and conjugacy is broken
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dgamma(1,1)
      mui[i] <- mu[xi[i]]
      y[i] ~ dexp(3*mui[i])
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  m = nimbleModel(code, data = list(y = rexp(4, 4)),
                  inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## non-exchangeable prior for tilde nodes
  code = nimbleCode({
    for(i in 1:4){
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], sd = 1)
      muTilde[i] ~ dnorm(mu0[i], sd = s0)
      mu0[i] ~ dnorm(0,1)
    }
    xi[1:4] ~ dCRP(1, 4)
    s0 ~ dhalfflat()
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4)))
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mcmc=buildMCMC(conf)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  
})  


test_that("Testing handling (including error detection) with non-standard CRP model specification",{
 
  n <- 20
  const <- list(n = n)
  inits <- list(xi = rep(1,n), muTilde = rnorm(n), conc = 1)
  data <- list(y = rnorm(n))
  tildeNames <- paste0("muTilde[", 1:n, "]")
  target <- paste0("xi[1:", n, "]")
  
  ## basic model to check results of findClusterNodes()
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n){
      muTilde[i] ~ dnorm(0,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## more complicated indexing to check results of findClusterNodes()
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i], 2]
    }
    for(i in 1:n){muTilde[i, 2] ~ dnorm(0,1)}
    
  })
  inits2 <- inits
  inits2$muTilde <- matrix(rnorm(n*2), n)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, ", 2]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(2, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## fewer tildeNodes than observations
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:(n-2)){muTilde[i] ~ dnorm(0,1)}
    
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "less than the number of potential clusters")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n-2, clusterNodeInfo$nTilde)
  
  ## indirect indexing; we don't have a good way to handle this.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[b[i]]
    }
    for(j in 1:n)
      b[j] <- xi[j]
    for(i in 1:n)
      muTilde[i] ~ dnorm(0,1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Detected that the CRP variable is used in some way not as an index")
  
  ## cluster ID as second index
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(muTilde[2, xi[i]], var = 1)
    }
    for(i in 1:n)
    {muTilde[2, i] ~ dnorm(0,1)}
  })
  inits2 <- inits
  inits2$muTilde <- rbind(rnorm(n), rnorm(n))
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[2, ", 1:n, "]"))
  expect_equal(2, clusterNodeInfo$numIndexes)
  expect_equal(2, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  
  ## clusterID as second index, additional nodes that are not clusterNodes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(muTilde[2, xi[i]], var = 1)
    }
    for(j in 1:2)
      for(i in 1:n)
      {muTilde[j, i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[2, ", 1:n, "]"))
  expect_equal(2, clusterNodeInfo$numIndexes)
  expect_equal(2, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  
  ## cluster ID in a function
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+1]
    }
    for(i in 1:(n+1))
    {muTilde[i] ~ dnorm(0,1)}
  })
  inits2$muTilde <- rnorm(n+1)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n+1), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## cluster ID in a function, no extraneous muTilde
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+1]
    }
    for(i in 2:(n+1))
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n+1), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## reordering of muTildes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[n-xi[i]+1]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", n:1, "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## function in index (not allowed)
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[n-i+1]]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  expect_error(conf <- configureMCMC(m), "findClusterNodes: Detected that a cluster parameter is indexed by a function")
  
  ## clusterNodes indexing doesn't begin at 1
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+2]
    }
    for(i in 3:(n+2))
    {muTilde[i] ~ dnorm(0,1)}
  })
  inits2$muTilde <- rnorm(n+2)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 3:(n+2), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## Extra nodes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:(2*n))
    {muTilde[i] ~ dnorm(0,1)}
  })
  inits2$muTilde <- rnorm(2*n)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## missing first cluster node: we no longer detect this situation because difficult to do with clusterNodes with more than one index.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 2:(n-2))
      muTilde[i] ~ dnorm(0,1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_warning(mcmc <- buildMCMC(conf), "sampler_CRP: The number of clusters based on the cluster parameters is less")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n-2), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n-3, clusterNodeInfo$nTilde)
  
  ## cluster node indexing shifted
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]+1]
    }
    for(i in 2:(n-2))
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "sampler_CRP: The number of clusters based on the cluster parameters is less")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 2:(n-2), "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(FALSE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n-3, clusterNodeInfo$nTilde)
  
  ## extra dependency on cluster nodes
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0, 1)}
    z ~ dnorm(muTilde[1], 1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Only the variables being clustered")

  
  ## cluster params depend on membership variable
  code=nimbleCode({
    for(i in 1:10) {
      muTilde[i] ~ dnorm(log(xi[1]), 1)  
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], 1)
    }
    xi[1:10] ~ dCRP(1 , size=10)
  })
  Inits=list(xi=rep(1, 10), muTilde=rep(0,10))
  Data=list(y=rnorm(10,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(mConf),
               'sampler_CRP: Detected that the CRP variable is used in some way not as an index')

       
  ## non related variable depends on cluster variable and membership variable
  code=nimbleCode({
    for(i in 1:10) {
      muTilde[i] ~ dnorm(0, 1)  
      mu[i] <- muTilde[xi[i]]
      y[i] ~ dnorm(mu[i], 1)
    }
    xi[1:10] ~ dCRP(1 , size=10)
    tau ~ dnorm(muTilde[xi[1]], 1)
  })
  Inits=list(xi=rep(1, 10), muTilde=rep(0,10), tau=1)
  Data=list(y=rnorm(10,0, 1))
  m <- nimbleModel(code, data=Data, inits=Inits)
  mConf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(mConf),
               'sampler_CRP: Detected unusual indexing')

  ## This is ok because y and x are same length.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i] ~ dnorm(0, var = s2[i]/lambda)
      s2[i] ~ dinvgamma(2, 1)
      y[i] ~ dnorm(mu[xi[i]],  var = s2[xi[i]])
      x[i] ~ dnorm(mu[xi[i]], 1)
    }
    lambda ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)
  })
  m <- nimbleModel(code, data=c(data, list(x = rnorm(n))), inits=inits, constants = const)
  mConf <- configureMCMC(m)
  mMCMC <- buildMCMC(mConf)

  ## Not ok because y and x are not the same length.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(alpha, n)
    for(i in 1:n){
      mu[i] ~ dnorm(0, var = s2[i]/lambda)
      s2[i] ~ dinvgamma(2, 1)
      y[i] ~ dnorm(mu[xi[i]],  var = s2[xi[i]])
    }
    for(i in 1:5) 
      x[i] ~ dnorm(mu[xi[i]], 1)
    lambda ~ dgamma(1, 1)
    alpha ~ dgamma(1, 1)
  })
  
  m <- nimbleModel(code, data=c(data, list(x = rnorm(5))), inits=inits, constants = const)
  mConf <- configureMCMC(m)
  expect_error(mMCMC <- buildMCMC(mConf), "sampler_CRP: Inconsistent indexing")
    
  ## Extraneous node that is ok.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n)
      muTilde[i] ~ dnorm(0,1)
    z ~ dnorm(muTilde[n+1], 1)
  })
  inits2$muTilde <- rnorm(n+1)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## Awkward trapping of observations with different distributions.
  ## This fails because we want the indexing of cluster parameters to be one contiguous block.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) 
      y[i] ~ dnorm(mu[i], var = 1)
    for(i in 1:(n-2))
        mu[i] <- muTilde[xi[i]]
    for(j in (n-1):n)
        mu[j] <- exp(muTilde[xi[j]])
    for(i in 1:n)
      muTilde[i] ~ dnorm(0, 1)
  })
  constSave <- const
  const$n <- 4
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "differing number of clusters indicated by")
  const <- constSave
  
  ## conjugate but observations not IID  
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2[i])
      mu[i] <- muTilde[xi[i]]
      s2[i] ~ dgamma(1,1)
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## conjugacy not detected because observations have multiple declarations
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:(n/2)) {
      y[i] ~ dnorm(mu[i], var = 1)
    }
    for(i in ((n/2)+1):n)
    {y[i] ~ dnorm(mu[i], var = 1)}
    for(i in 1:n)
    {mu[i] <- muTilde[xi[i]]}
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## observations not independent
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    y[1] ~ dnorm(mu[1], var = s2[1])
    for(i in 2:n)
      y[i] ~ dnorm(mu[i]+y[i-1], var = s2[i])
    for(i in 1:n) {
      mu[i] <- muTilde[xi[i]]
      s2[i] ~ dgamma(1,1)
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Variables being clustered must be conditionally independent.")
  
  ## cluster nodes not independent
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    muTilde[1] ~ dnorm(0, 1)
    for(i in 2:n)
    {muTilde[i] ~ dnorm(muTilde[i-1],1)}
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: cluster parameters must be independent across clusters")
  
  ## cluster nodes not exchangeable so non-conjugate
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:(n-1) )
    {muTilde[i] ~ dnorm(mu0[i],1)}
    muTilde[n] ~ dgamma(1,1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("muTilde[", 1:n, "]"))
  expect_equal(TRUE, clusterNodeInfo$targetIsIndex)
  expect_equal(FALSE, clusterNodeInfo$targetIndexedByFunction)
  expect_equal(1, clusterNodeInfo$numIndexes)
  expect_equal(1, clusterNodeInfo$indexPosition)
  expect_equal(n, clusterNodeInfo$nTilde)
  
  ## cluster membership variables not independent of cluster parameters
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc + muTilde[1], n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n)
    {muTilde[i] ~ dnorm(0,1)}
    
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Only the variables being clustered can depend")
  
  ## cluster membership variables not independent of cluster parameters
  ## This is not detected by the check for independence but by use of 'xi' as non-index
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = 1)
      mu[i] <- muTilde[xi[i]]
      tmp[i] ~ dnorm(0,1)
    }
    for(i in 1:n)
      muTilde[i] ~ dnorm(tmp[xi[i]],1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "Only the variables being clustered can depend")
  
  inits$s2Tilde <- rep(1, n)
  
  ## Conjugate normal-invgamma
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0, var = s2Tilde[i])
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:n, "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", 1:n, "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(rep(TRUE, 2), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(rep(n,2), clusterNodeInfo$nTilde)
  
  ## nTilde < n for one of the cluster parameters. Note confusing error message.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i]]
    }
    kappa ~ dgamma(1,1)
    for(i in 1:n) 
      muTilde[i] ~ dnorm(0, var = s2Tilde[i]/kappa)
    for(i in 1:(n-1))
      s2Tilde[i] ~ dinvgamma(1,1)
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: In a model with multiple cluster parameters, the number")
  
  ## nTilde < n 
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i]]
    }
    kappa ~ dgamma(1,1)
    for(i in 1:(n-1)) {
      muTilde[i] ~ dnorm(0,var = s2Tilde[i]/kappa)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "less than the number of potential clusters")
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:(n-1), "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", 1:(n-1), "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(rep(TRUE, 2), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(c(n-1, n-1), clusterNodeInfo$nTilde)
  
  ## CRP variable used in multiple indices; disallowing this.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[xi[i]])
      mu[i] <- muTilde[xi[i],xi[i]]
    }
    for(i in 1:n) {
      for(j in 1:n)
        {muTilde[i,j] ~ dnorm(0,var=s2Tilde[i]/3)}
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  inits2 <- inits
  inits2$muTilde <- matrix(rnorm(n^2),n)
  m <- nimbleModel(code, data = data, constants = const, inits = inits2)
  expect_error(conf <- configureMCMC(m), "CRP variable used multiple times")
  
  ## weird ordering of muTilde/s2Tilde but should be ok
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[n-xi[i]+1])
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,var=s2Tilde[n-i+1]/3)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_invgamma_dnorm")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:n, "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", n:1, "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(c(FALSE, TRUE), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(rep(n,2), clusterNodeInfo$nTilde)
  
  ## s2Tildes in different order than muTildes so not conjugate.
  ## CRP_sampler would INCORRECT for this because can't sample from distr of an s2Tilde given the muTilde that depends on it.
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(mu[i], var = s2Tilde[n-xi[i]+1])
      mu[i] <- muTilde[xi[i]]
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,var=s2Tilde[i]/3)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: cluster parameters must be independent across clusters")
    
  ## Non-conjugate, bivariate
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      y[i] ~ dnorm(muTilde[xi[i]], var = exp(s2Tilde[xi[i]]))
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,1)
      s2Tilde[i] ~ dgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  conf <- configureMCMC(m)
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  clusterNodeInfo <- nimble:::findClusterNodes(m, target)
  expect_equal(clusterNodeInfo$clusterNodes[[2]], paste0("muTilde[", 1:n, "]"))
  expect_equal(clusterNodeInfo$clusterNodes[[1]], paste0("s2Tilde[", 1:n, "]"))
  expect_equal(c(1,1), clusterNodeInfo$numIndexes)
  expect_equal(c(1,1), clusterNodeInfo$indexPosition)
  expect_equal(rep(TRUE, 2), clusterNodeInfo$targetIsIndex)
  expect_equal(rep(FALSE, 2), clusterNodeInfo$targetIndexedByFunction)
  expect_equal(rep(n,2), clusterNodeInfo$nTilde)
  
  ## cross clustering
  
  data$y <- matrix(rnorm(n^2), n)
  
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      for(j in 1:n)
        y[i,j] ~ dnorm(muTilde[xi[i]], var = s2Tilde[xi[j]])
    }
    for(i in 1:n) {
      muTilde[i] ~ dnorm(0,1)
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  expect_error(conf <- configureMCMC(m), "findClusterNodes: found cluster membership parameters that use different indexing variables")
  
  inits$muTilde <- matrix(rnorm(n^2), n)
  code <- nimbleCode({
    xi[1:n] ~ dCRP(conc, n)
    for(i in 1:n) {
      for(j in 1:n)
        {y[i,j] ~ dnorm(muTilde[xi[i],xi[j]], var = s2Tilde[xi[i]])}
    }
    for(i in 1:n)  {
      for(j in 1:n)
        {muTilde[i,j] ~ dnorm(0,1)}
      s2Tilde[i] ~ dinvgamma(1,1)
    }
  })
  m <- nimbleModel(code, data = data, constants = const, inits = inits)
  expect_error(conf <- configureMCMC(m), "CRP variable used multiple times in")

  ## Various additional cases based on more general BNP functionality

  ## Basic more general case
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })
  
  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  expect_identical(cn$clusterNodes[[1]], c(matrix(model$expandNodeNames('thetaTilde'), J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))

  ## Basic case of truncated parameters
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
          }}
      for(i in 1:n2) {
          for(j in 1:J) {
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  n2 <- 4
  J <- 3
  constants <- list(n = n, n2 = n2, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n2), n2, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  expect_identical(cn$clusterNodes[[1]], c(matrix(model$expandNodeNames('thetaTilde'), J, n2, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "is less than the number of potential")
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n2*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n2, each = J))

  ## Truncated and intermediate
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(theta[i, j], 1)
              theta[i, j] <- thetaTilde[xi[i], j]
          }}
      for(i in 1:n2)
          for(j in 1:J)
              thetaTilde[i, j] ~ dnorm(0, 1)
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  n2 <- 4
  J <- 3
  constants <- list(n = n, n2 = n2, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  nodes <- nodes[nodes %in% model$getNodeNames()]
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n2, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_warning(mcmc <- buildMCMC(conf), "is less than the number of potential")
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n2*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n2, each = J))

  ## Column instead of row
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[j, i] ~ dnorm(thetaTilde[j, xi[i]], 1)
              thetaTilde[j, i] ~ dnorm(0, 1)
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),J,n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), J, n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  expect_identical(cn$clusterNodes[[1]], model$expandNodeNames('thetaTilde'))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))

  ## Column in instead of row, different loop ordering
  code <- nimbleCode({
      for(j in 1:J) {
          for(i in 1:n) {
              y[j, i] ~ dnorm(thetaTilde[j, xi[i]], 1)
              thetaTilde[j, i] ~ dnorm(0, 1)
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),J,n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), J, n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  expect_identical(cn$clusterNodes[[1]], model$expandNodeNames('thetaTilde'))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))

  ## Column instead of row, intermediate nodes, indexes swapped
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[j, i] ~ dnorm(theta[j, i], 1)
              theta[j, i] <- thetaTilde[xi[i], j]
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),J,n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))


  ## index offset
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i]+1, j], 1)
          }}
      for(i in 2:(n+1)) {
          for(j in 1:J) {
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*(n+1)), n+1, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n+1, byrow = TRUE)[,-1]))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))


  ## index offset, with intermediate
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(theta[i, j], 1)
              theta[i,j] <- thetaTilde[xi[i]+1, j]
          }}
      for(i in 2:(n+1)) {
          for(j in 1:J) {
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*(n+1)), n+1, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n+1, byrow = TRUE)[,-1]))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))

  ## multiple obs, but only single thetaTilde in each group
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i]], 1)
          }
          thetaTilde[i] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  expect_identical(cn$clusterNodes[[1]], model$expandNodeNames('thetaTilde'))
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)
  
  ## Detection of differing number of cluster parameters
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[i]], s2tilde[xi[i]])
          s2tilde[i] ~ dunif(0, 10)
      }
      for(i in 1:n2) 
          thetaTilde[i] ~ dnorm(0, 1)
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  n2 <- 4
  J <- 3
  constants <- list(n = n, n2 = n2, J = J)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n2), s2tilde = runif(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: In a model with multiple cluster parameters")

  ## Detection of differing number of cluster parameters, complicated indexing
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[i]], s2tilde[xi[i]+1])
          s2tilde[i] ~ dunif(0, 10)
          thetaTilde[i] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n), s2tilde = runif(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: In a model with multiple cluster parameters")

  ## Another variation where we have intermediate nodes and dependency in cluster nodes
  ## clustering of deterministic nodes - not allowed
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(theta[i], sd = sigma[i])
          theta[i] <- thetaTilde[xi[i]]
          sigma[i] <- sigmaTilde[xi[i]]
          sigmaTilde[i] <- 1 / tauTilde[i]
          thetaTilde[i] ~ dnorm(mu, var = sigmaTilde[i])
          tauTilde[i] ~ dgamma(a, b)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })
  n <- 5
  J <- 2
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n), tauTilde = rgamma(n, 1, 1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "findClusterNodes: detected that deterministic nodes are being clustered")

  ## Model A: dnorm obs, dmnorm prior
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J)
              y[i,j] ~ dnorm(thetaTilde[xi[i],j], 1)
          thetaTilde[i, 1:J] ~ dmnorm(mn[1:J], iden[1:J,1:J])
      }
      mn[1:J] <- mu*ones[1:J]
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J),
                iden = diag(3), mu = 0)
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)

  ## Model B: separate prior specifications
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:2) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
          }
          thetaTilde[i, 1] ~ dnorm(0, 1)
          thetaTilde[i, 2] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n, J=J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow =  TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))


  ## Model C: dmnorm obs, dmnorm prior
  code <- nimbleCode({
      for(i in 1:n) {
          y[i, 1:J] ~ dmnorm(thetaTilde[xi[i],1:J], iden[1:J,1:J])
          thetaTilde[i, 1:J] ~ dmnorm(mn[1:J], iden[1:J,1:J])
      }
      mn[1:J] <- mu*ones[1:J]
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J),
                iden = diag(3), mu = 0)
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dmnorm_dmnorm")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)


  ## Model D: dmnorm obs, dnorm prior
  code <- nimbleCode({
      for(i in 1:n) {
          y[i, 1:J] ~ dmnorm(thetaTilde[xi[i],1:J], iden[1:J,1:J])
          for(j in 1:J)
              thetaTilde[i, j] ~ dnorm(0,1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J),
                iden = diag(J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))

  ## Model E: dmnorm obs, separate prior declarations
  code <- nimbleCode({
      for(i in 1:n) {
          y[i, 1:2] ~ dmnorm(thetaTilde[xi[i],1:2], iden[1:2,1:2])
          thetaTilde[i, 1] ~ dnorm(0,1)
          thetaTilde[i, 2] ~ dnorm(0,1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J),
                iden = diag(J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))

  ## Model F: dnorm obs, separate priors
  code <- nimbleCode({
      for(i in 1:n) {
          y1[i] ~ dnorm(thetaTilde1[xi[i]], 1)
          y2[i] ~ dnorm(thetaTilde2[xi[i]], 1)
          thetaTilde1[i] ~ dnorm(mu, sigma)
          thetaTilde2[i] ~ dnorm(mu, sigma)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y1 = rnorm(n), y2 = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde1 = rnorm(n), thetaTilde2 = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde1')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde1')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)

  ## Model F: non conjugate, non-identical priors
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
      }
      thetaTilde[1] ~ dnorm(mu, sigma)
      for(i in 2:n) 
          thetaTilde[i] ~ dgamma(1,1)
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), as.integer(c(2:n, 1)))

  ## Model F: non conjugate, non-identical data
  ## doesn't detect presence of some conjugacies because only check first set
  ## weird setup causes us to say we can't handle this case
  code <- nimbleCode({
      for(i in 3:n) {
          y[i] ~ dt(thetaTilde[xi[i]], 1, 1)
      }
      for(i in 1:2)
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
      for(i in 1:n) 
          thetaTilde[i] ~ dnorm(0,1)
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: In a model with multiple cluster parameters")

  ## Model F: non conjugate, non-identical data nodes
  code <- nimbleCode({
      for(i in 2:n) {
          y[i] ~ dt(thetaTilde[xi[i]], 1, 1)
      }
      y[1] ~ dnorm(thetaTilde[xi[1]], 1)
      for(i in 1:n) 
          thetaTilde[i] ~ dnorm(0,1)
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Detected unusual indexing")

  ## another case of mixed distributions
  code <- nimbleCode({
      for(i in 1:(n-1)) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      for(j in 1:J) {
          y[n, j] ~ dt(thetaTilde[xi[n], j], 1 ,1)
          thetaTilde[n, j] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Detected unusual indexing")

  ## Model G: separate declarations for obs
  code <- nimbleCode({
      for(i in 1:n) {
          y1[i] ~ dnorm(thetaTilde[xi[i], 1], 1)
          y2[i] ~ dnorm(thetaTilde[xi[i], 2], 1)
          for(j in 1:2)
              thetaTilde[i, j] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n)
  data <- list(y1 = rnorm(n), y2 = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(n*J),n,J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes[1:n])
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))


  ## Model H: separate obs and prior declarations
  code <- nimbleCode({
      for(i in 1:n) {
          y1[i] ~ dnorm(thetaTilde[xi[i], 1], 1)
          y2[i] ~ dnorm(thetaTilde[xi[i], 2], 1)
          thetaTilde[i, 1] ~ dnorm(0, 1)
          thetaTilde[i, 2] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n)
  data <- list(y1 = rnorm(n), y2 = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(n*J),n,J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes[1:n])
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))

  ## Model I: separate declarations for obs, dmnorm for prior
  code <- nimbleCode({
      for(i in 1:n) {
          y1[i] ~ dnorm(thetaTilde[xi[i], 1], 1)
          y2[i] ~ dnorm(thetaTilde[xi[i], 2], 1)
          thetaTilde[i, 1:2] ~ dmnorm(mn[1:2], iden[1:2,1:2])
      }
      mn[1:2] <- mu*ones[1:2]
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n)
  data <- list(y1 = rnorm(n), y2 = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n), iden = diag(1, 2),
                thetaTilde = matrix(rnorm(n*J), n, J), ones = rep(1,2))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)

  ## model J: mixed indexing of two cluster vars
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i,j] ~ dnorm(theta[i], sd = sigma[i,j])        
              sigma[i,j] <- sigmaTilde[xi[i], j]
              sigmaTilde[i,j] ~ dinvgamma(a, b)
          }
          theta[i] <- thetaTilde[xi[i]]
          thetaTilde[i] ~ dnorm(mu, phi)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J), n, J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n), sigmaTilde = matrix(rgamma(n*J, 1, 1), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], nodes)
  nodes <- model$expandNodeNames('sigmaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, c(2L, 1L))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J+1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('thetaTilde'), conf$getSamplers('sigmaTilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*(J+1)))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), c(1:n, rep(1:n, each = J)))


  ## Model K: more complicated intermediate nodes
  ## We do not have normal-gamma conj set up; this is ok.
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(theta[i], sd = sigma[i])
          theta[i] <- thetaTilde[xi[i]]
          sigma[i] <- 1 / tau[i]
          tau[i] <- tauTilde[xi[i]]
          thetaTilde[i] ~ dnorm(mu, var = sigmaTilde[i])
          sigmaTilde[i] <- 1 / tauTilde[i]
          tauTilde[i] ~ dgamma(a, b)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })
  n <- 5
  J <- 2
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n), tauTilde = rgamma(n, 1, 1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], nodes)
  nodes <- model$expandNodeNames('tauTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(3))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('thetaTilde'), conf$getSamplers('tauTilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))

  ## 3-d case with non-variable 3rd index
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i,j] ~ dnorm(thetaTilde[xi[i], 2, j] , var = 1)
          }
      }
      for(i in 1:n) {
          for(j in 1:J) {
              for(k in 1:2) {
                  thetaTilde[i, k, j] ~ dnorm(0,1)
              }
          }
      }
      xi[1:n] ~ dCRP(1, size=n)
  })

  n <- 5
  J <- 3
  constants  <- list(n = n, J = J) 
  inits <- list(xi = rep(1, n), 
                thetaTilde = array(0, c(n,2,J)))
  y <- matrix(0, nrow = n , ncol= J) 
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(t(array(nodes, c(n, 2, J))[,2,])))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')[16:30]
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))


  ## 3-d case with full third index
  code <- nimbleCode({
      for(i in 1:n) {
          for(k in 1:2) {
              for(j in 1:J) {
                  y[k,i,j] ~ dnorm(thetaTilde[k, xi[i], j] , var = 1)
              }
          }
      }
      for(i in 1:n) {
          for(j in 1:J) {
              for(k in 1:2) {
                  thetaTilde[k, i, j] ~ dnorm(0,1)
              }
          }
      }
      xi[1:n] ~ dCRP(1, size=n)
  })

  n <- 5
  J <- 3
  constants  <- list(n = n, J = J) 
  inits <- list(xi = rep(1, n), 
                thetaTilde = array(0, c(2, n, J)))
  y <- array(0, c(2, n, J))
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- c(matrix(model$expandNodeNames('thetaTilde[1:2, 1, 1:3]'), J, 2, byrow = TRUE),
             matrix(model$expandNodeNames('thetaTilde[1:2, 2, 1:3]'), J, 2, byrow = TRUE),
             matrix(model$expandNodeNames('thetaTilde[1:2, 3, 1:3]'), J, 2, byrow = TRUE),
             matrix(model$expandNodeNames('thetaTilde[1:2, 4, 1:3]'), J, 2, byrow = TRUE),
             matrix(model$expandNodeNames('thetaTilde[1:2, 5, 1:3]'), J, 2, byrow = TRUE))
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, 6L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, 2*J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J*2))

  ## index is a function - not allowed
  code <- nimbleCode({
      for(i in 1:3) {
          for(j in 1:4) {
              y[i,j] ~ dnorm(thetaTilde[xi[3-i+1], j] , var = 1) 
              thetaTilde[i, j] ~ dnorm(0,1)
          }
      }
      xi[1:3] ~ dCRP(1, size=3)
  })

  inits <- list(xi = c(1, 1, 1), 
                thetaTilde = matrix(0, nrow=3,  ncol=4))
  y <- matrix(5, nrow=3, ncol=4) 
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
               "findClusterNodes: Detected that a cluster parameter is indexed by a function")

  ## index is a function - not allowed
  code <- nimbleCode({
      for(i in 1:3) {
          for(j in 1:4) {
              y[i,j] ~ dnorm(thetaTilde[xi[i+j], j] , var = 1) 
              thetaTilde[i, j] ~ dnorm(0,1)
          }
      }
      xi[1:7] ~ dCRP(1, size=7)
  })

  inits <- list(xi = rep(1,7), 
                thetaTilde = matrix(0, nrow=3,  ncol=4))

  y <- matrix(5, nrow=3, ncol=4) 
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
               "findClusterNodes: Detected that a cluster parameter is indexed by a function")

  ## Cases of multiple use of CRP node 

  ## use of xi[i] and xi[j] in same statement - not allowed
  code <- nimbleCode({
      for(i in 1:3) {
          for(j in 1:4) {
              y[i,j] ~ dnorm(thetaTilde[xi[i], xi[j]] , var = 1) 
              thetaTilde[i, j] ~ dnorm(0,1)
          }
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  inits <- list(xi = rep(1,7), 
                thetaTilde = matrix(0, nrow=3,  ncol=4))

  y <- matrix(5, nrow=3, ncol=4) 
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
               "findClusterNodes: CRP variable used multiple times")

  ## use of xi[i] and xi[i] in same parameter - not allowed
  code <- nimbleCode({
      for(i in 1:4) {
          for(j in 1:4) {
              y[i,j] ~ dnorm(thetaTilde[xi[i], xi[i]] , var = 1) 
              thetaTilde[i, j] ~ dnorm(0,1)
          }
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  inits <- list(xi = rep(1,7), 
                thetaTilde = matrix(0, nrow=4,  ncol=4))

  y <- matrix(5, nrow=4, ncol=4)
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
               "findClusterNodes: CRP variable used multiple times")

  ## xi[i] and xi[j] in separate parameters - not allowed
  code <- nimbleCode({
      for(i in 1:4) {
          for(j in 1:4) {
              y[i,j] ~ dnorm(thetaTilde[xi[i]], var = s2Tilde[xi[j]])
          }
      }
      for(i in 1:4)
          thetaTilde[i] ~ dnorm(0,1)
      for(i in 1:4)
          s2Tilde[i] ~ dnorm(0,1)
      xi[1:4] ~ dCRP(1, size=4)
  })
  inits <- list(xi = rep(1,4))
  y <- matrix(5, nrow=4, ncol=4)
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
               "findClusterNodes: found cluster membership parameters that use different indexing variables")

  ## xi[i] used twice in same parameter legitimately
  ## conjugate but we are not set up to use this conjugacy
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(b0[xi[i]] + b1[xi[i]]*x[i], var = 1)
      }
      for(i in 1:n) {
          b0[i] ~ dnorm(0,1)
          b1[i] ~ dnorm(0,1)
      }
      xi[1:n] ~ dCRP(1, size=n)
  })
  n  <- 5
  constants <- list(n = n)
  data = list(y = rnorm(n))
  inits = list(x = rnorm(n), xi = rep(1,n))
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('b0')
  expect_identical(cn$clusterNodes[[1]], nodes)
  nodes <- model$expandNodeNames('b1')
  expect_identical(cn$clusterNodes[[2]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('b0'), conf$getSamplers('b1'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))

  ## xi[i] used twice in same parameter, variation
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(beta[xi[i], 1] + beta[xi[i], 2]*x[i], var = 1)
      }
      for(i in 1:n) {
          beta[i,1] ~ dnorm(0,1)
          beta[i,2] ~ dnorm(0,1)
      }
      xi[1:n] ~ dCRP(1, size = n)
  })
  n  <- 5
  constants <- list(n = n)
  data = list(y = rnorm(n))
  inits = list(x = rnorm(n), xi = rep(1,n))
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('beta')
  expect_identical(cn$clusterNodes[[1]], nodes[1:n])
  expect_identical(cn$clusterNodes[[2]], nodes[(n+1):(2*n)])
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('beta')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))

  ## xi[i] used twice in same parameter, another variation
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(beta[xi[i], 1] + beta[xi[i], 2]*x[i], var = 1)
      }
      for(i in 1:n) 
          for(j in 1:2) 
              beta[i,j] ~ dnorm(0,1)
      xi[1:n] ~ dCRP(1, size=n)
  })
  n  <- 5
  constants <- list(n = n)
  data = list(y = rnorm(n))
  inits = list(x = rnorm(n), xi = rep(1,n))
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('beta')
  expect_identical(cn$clusterNodes[[1]], nodes[1:n])
  expect_identical(cn$clusterNodes[[2]], nodes[(n+1):(2*n)])
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('beta')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))


  ## use of inprod
  ## conj not detected as we are only set up for b0 + x[i]*b1[xi[i]] type conjugacy
  ## not b0[xi[i]]+x[i]*b1[xi[i]]
  ## need to think more about when the nonidentity dnorm-dnorm conjugacy would be used
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(inprod(beta[1:J, xi[i]], x[i,1:J]), var = 1)
      }
      for(i in 1:n)
          for(j in 1:J)
              beta[j, i] ~ dnorm(0,1)
      xi[1:n] ~ dCRP(1, size = n)
  })
  n  <- 5
  J <- 2
  constants <- list(n = n, J = J)
  data = list(y = rnorm(n))
  inits = list(x = matrix(rnorm(n*J),n,J), xi = rep(1,n))
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('beta')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('beta')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))


  ## a variation on the previous case, changing the prior
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(inprod(beta[1:J, xi[i]], x[i,1:J]), var = 1)
      }
      for(i in 1:n)
          beta[1:J, i] ~ dmnorm(z[1:J], pr[1:J,1:J])
      xi[1:n] ~ dCRP(1, size=n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n, J = J)
  data = list(y = rnorm(n))
  inits = list(x = matrix(rnorm(n*J),n,J), xi = rep(1,n), pr = diag(J))
  model <- nimbleModel(code, data = data, inits = inits, constants = constants)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('beta')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, as.integer(1))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('beta')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)

  ## non-independent observations within a group - not allowed
  code <- nimbleCode({
      for(i in 1:4) {
          for(j in 2:3) {
              y[i,j] ~ dnorm(mu[xi[i]] + y[i,j-1], 1)
          }
      }
      for(i in 1:4) {
          mu[i] ~ dnorm(0,1)
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  data = list(y =matrix( rnorm(4*3), 4 ,3))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Variables being clustered must be conditionally independent")

  ## another non-indep of observations within a group - not allowed
  code <- nimbleCode({
      for(i in 1:4) {
          for(j in 2:3) {
              y[i,j] ~ dnorm(mu[xi[i]], exp(y[i,j-1]))
          }
      }
      for(i in 1:4) {
          mu[i] ~ dnorm(0,1)
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  data = list(y =matrix( rnorm(4*3), 4 ,3))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Variables being clustered must be conditionally independent")

  ## observations dependent across group - not allowed
  code <- nimbleCode({
      for(i in 2:4) {
          for(j in 1:3) {
              y[i,j] ~ dnorm(mu[xi[i]] + y[i-1,j], 1)
          }
      }
      for(i in 1:4) {
          mu[i] ~ dnorm(0,1)
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  data = list(y =matrix( rnorm(4*3), 4 ,3))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Variables being clustered must be conditionally independent")

  ## observations dependent across group, second case - not allowed
  code <- nimbleCode({
      for(i in 1:3) {
          for(j in 1:3) {
              y[i,j] ~ dnorm(mu[xi[i]] + y[i+1,j], 1)
          }
      }
      for(i in 1:4) {
          mu[i] ~ dnorm(0,1)
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  data = list(y =matrix( rnorm(4*3), 4 ,3))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  ## errors with lack of independence but error doesn't distinguish within vs. between group
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Variables being clustered must be conditionally independent")

  ## non-independence across groups, third case - not allowed
  code <- nimbleCode({
      for(i in 1:4) {
          y[i,1] ~ dnorm(mu[xi[i]], 1)
          for(j in 2:3) {
              y[i,j] ~ dnorm(mu[xi[i]] + y[i,j-1], 1)
          }
      }
      for(i in 1:4) {
          mu[i] ~ dnorm(0,1)
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  data = list(y =matrix( rnorm(4*3), 4 ,3))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  ## errors with lack of independence but error doesn't distinguish within vs. between group
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Variables being clustered must be conditionally independent")

  ## y dependence in multivariate way
  code <- nimbleCode({
      for(i in 1:n) {
          y[i,1:2] ~ dmnorm(mu[xi[i], 1:2], cov = sigma[xi[i],1:2,1:2])
          mu[i, 1:2] ~ dmnorm(mu0[1:2], iden[1:2,1:2])
          sigma[i,1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
      }
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  data = list(y =matrix( rnorm(n*2), n ,2))
  sigma <- array(0, c(n, 2, 2))
  for(i in 1:n)
      sigma[i, 1:2, 1:2] <- diag(2)
  inits = list(xi = rep(1,n), iden = diag(2), sigma = sigma, S = diag(2))
  model <- nimbleModel(code, constants = list(n = n), data = data, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], nodes)
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))

  ## various dmnorm-invwish cases
  sigma <- array(0, c(n, 2, 2))
  for(i in 1:n)
      sigma[i, 1:2, 1:2] <- diag(2)

  ## single obs per clusterID
  code <- nimbleCode({
      for(i in 1:n) {
          y[i,1:2] ~ dmnorm(mu[xi[i], 1:2], cov = sigma[xi[i],1:2,1:2])
          mu[i, 1:2] ~ dmnorm(mu0[1:2], cov = sigmaAux[i,1:2,1:2])
          sigmaAux[i,1:2,1:2] <- sigma[i,1:2,1:2]/kappa
          sigma[i,1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
      }
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  data = list(y =matrix( rnorm(n*2), n ,2))
  inits = list(xi = rep(1,n), S = diag(2), sigma = sigma, kappa = 1)
  model <- nimbleModel(code, constants = list(n = n), data = data, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], nodes)
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dmnorm_invwish_dmnorm")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))

  ## standard case with multiple obs
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:2) {
              y[i, j,1:2] ~ dmnorm(mu[xi[i], j, 1:2], cov = sigma[xi[i], j, 1:2,1:2])
              mu[i, j, 1:2] ~ dmnorm(mu0[1:2], cov = sigmaAux[i, j, 1:2,1:2])
              sigmaAux[i, j, 1:2,1:2] <- sigma[i, j, 1:2,1:2]/kappa
              sigma[i, j, 1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
          }}
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  sigma <- array(0, c(n, 2, 2, 2))
  for(i in 1:n)
      for(j in 1:2) 
      sigma[i, j, 1:2, 1:2] <- diag(2)
  data = list(y =array(rnorm(n*2*2), c(n, 2, 2)))
  inits = list(xi = rep(1,n), kappa = 1, sigma = sigma, S = diag(2))
  model <- nimbleModel(code, constants = list(n=n), data = data, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, 2, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, 2, n, byrow = TRUE)))
  
  expect_identical(cn$numNodesPerCluster, rep(2L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dmnorm_invwish_dmnorm")
  expect_identical(crpSampler$nObsPerClusID, 2)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(2))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = 2), 2))


  ## not IID across clusters
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:2) {
              y[i, j,1:2] ~ dmnorm(mu[xi[i], j, 1:2], cov = sigma[xi[i], j, 1:2,1:2])
              mu[i, j, 1:2] ~ dmnorm(mu0[i, 1:2], cov = sigmaAux[i, j, 1:2,1:2])
              sigmaAux[i, j, 1:2,1:2] <- sigma[i, j, 1:2,1:2]/kappa
              sigma[i, j, 1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
          }}
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  data = list(y =array(rnorm(n*2*2), c(n, 2, 2)))
  inits = list(xi = rep(1,n), kappa = 1, S = diag(2), sigma = sigma)
  model <- nimbleModel(code, data = data, inits = inits, constants = list(n = n))
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, 2, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, 2, n, byrow = TRUE)))
  
  expect_identical(cn$numNodesPerCluster, rep(2L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 2)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(2))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = 2), 2))


  ## IID across but not within clusters
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:2) {
              y[i, j,1:2] ~ dmnorm(mu[xi[i], j, 1:2], cov = sigma[xi[i], j, 1:2,1:2])
              mu[i, j, 1:2] ~ dmnorm(mu0[j, 1:2], cov = sigmaAux[i, j, 1:2,1:2])
              sigmaAux[i, j, 1:2,1:2] <- sigma[i, j, 1:2,1:2]/kappa
              sigma[i, j, 1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
          }}
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  data = list(y =array(rnorm(n*2*2), c(n, 2, 2)))
  inits = list(xi = rep(1,n), kappa = 1, sigma = sigma, S = diag(2))
  model <- nimbleModel(code, data = data, inits = inits, constants = list(n=n))

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, 2, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, 2, n, byrow = TRUE)))
  
  expect_identical(cn$numNodesPerCluster, rep(2L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dmnorm_invwish_dmnorm")
  expect_identical(crpSampler$nObsPerClusID, 2)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(2))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = 2), 2))


  ## no dependence on sigma in mu prior
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:2) {
              y[i, j,1:2] ~ dmnorm(mu[xi[i], j, 1:2], cov = sigma[xi[i], j, 1:2,1:2])
              mu[i, j, 1:2] ~ dmnorm(mu0[1:2], cov = pr[1:2,1:2])
              sigma[i, j, 1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
          }}
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  data = list(y =array(rnorm(n*2*2), c(n, 2, 2)))
  inits = list(xi = rep(1,n), sigma = sigma, pr = diag(2), S = diag(2))
  model <- nimbleModel(code, data = data, inits = inits, constants = list(n=n))

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, 2, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, 2, n, byrow = TRUE)))
  
  expect_identical(cn$numNodesPerCluster, rep(2L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 2)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(2))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = 2), 2))

  ## only one sigma per cluster
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:2) {
              y[i, j, 1:2] ~ dmnorm(mu[xi[i], j, 1:2], cov = sigma[xi[i], 1:2,1:2])
              mu[i, j, 1:2] ~ dmnorm(mu0[1:2], cov = sigmaAux[i, 1:2,1:2])
          }
          sigmaAux[i, 1:2,1:2] <- sigma[i, 1:2,1:2]/kappa
          sigma[i, 1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
      }
      xi[1:n] ~ dCRP(1, size=n)
  })
  sigma <- array(0, c(n, 2, 2))
  for(i in 1:n)
      sigma[i, 1:2, 1:2] <- diag(2)
  n <- 5
  data = list(y =array(rnorm(n*2*2), c(n, 2, 2)))
  inits = list(xi = rep(1,n), sigma = sigma, S = diag(2), kappa = 1)
  model <- nimbleModel(code, data = data, inits = inits, constants = list(n = n))

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, 2, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, c(1L, 2L))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 2)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*3))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), c(rep(1:n, each = 2), 1:n))

  ## only one mu and one sigma per cluster
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:2) 
              y[i, j,1:2] ~ dmnorm(mu[xi[i], 1:2], cov = sigma[xi[i], 1:2,1:2])
          mu[i, 1:2] ~ dmnorm(mu0[1:2], cov = sigmaAux[i, 1:2,1:2])
          sigmaAux[i, 1:2,1:2] <- sigma[i, 1:2,1:2]/kappa
          sigma[i, 1:2,1:2] ~ dinvwish(S[1:2,1:2], nu)
      }
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  data = list(y =array(rnorm(n*2*2), c(5, 2, 2)))
  inits = list(xi = rep(1,n), kappa = 1, S = diag(2), sigma = sigma)
  model <- nimbleModel(code, data = data, inits = inits, constants = list(n=n))
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[2]], nodes)
  nodes <- model$expandNodeNames('sigma')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 2)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- c(conf$getSamplers('mu'), conf$getSamplers('sigma'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))

  ## y dependence in multivariate way
  code <- nimbleCode({
      for(i in 1:n) {
          y[i,1:2] ~ dmnorm(mu[xi[i], 1:2], cov = sigma[1:2,1:2])
          mu[i, 1:2] ~ dmnorm(mu0[1:2], iden[1:2,1:2])
      }
      xi[1:n] ~ dCRP(1, size=n)
  })
  n <- 5
  data = list(y =matrix( rnorm(n*2), n ,2))
  inits = list(xi = rep(1,n), iden = diag(2), sigma = diag(2))
  model <- nimbleModel(code, data = data, inits = inits, constants = list(n=n))
  
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('mu')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dmnorm_dmnorm")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('mu')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)


  ## clusters not independent - not allowed
  code <- nimbleCode({
      for(i in 1:4) 
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
      thetaTilde[1] ~ dnorm(a,1)
      for(i in 2:4)
          thetaTilde[i] ~ dnorm(thetaTilde[i-1], 1)
      xi[1:4] ~ dCRP(1, size=4)
      a ~ dunif(0,1)
  })
  data = list(y = rnorm(4))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: cluster parameters must be independent across clusters")

  ## clusters not indep - alternative case - not allowed
  code <- nimbleCode({
      for(i in 1:4) 
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
      thetaTilde[4] ~ dnorm(0,1)
      for(i in 1:3)
          thetaTilde[i] ~ dnorm(thetaTilde[i+1], 1)
      xi[1:4] ~ dCRP(1, size=4)
  })
  data = list(y = rnorm(4))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: cluster parameters must be independent across clusters")


  ## clusters not indep, with mv declaration - not allowed
  ## errors when trying to wrap sampler because there is only
  ## one cluster node. Might want have configureMCMC
  ## give up on wrapping and then have error occur in buildMCMC().
  code <- nimbleCode({
      for(i in 1:4) 
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
      thetaTilde[1:4] ~ dmnorm(z[1:4], iden[1:4,1:4]) # formally indep but not known to us
      xi[1:4] ~ dCRP(1, size=4)
  })
  data = list(y = rnorm(4))
  inits = list(xi = rep(1,4), iden = diag(4))
  model <- nimbleModel(code, data = data, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
               "Cannot determine wrapped sampler for cluster parameter")

  ## clusters indep G0 but not IID
  code <- nimbleCode({
      for(i in 1:4)  {
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
          thetaTilde[i] ~ dnorm(i, 1)
      }
      xi[1:4] ~ dCRP(1, size=4)
  })
  n <- 4
  data = list(y = rnorm(4))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:4]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)


  ## clusters indep G0 but not IID, case 2
  code <- nimbleCode({
      for(i in 1:4) 
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
      for(i in 1:3) 
          thetaTilde[i] ~ dnorm(0, 1)
      thetaTilde[4] ~ dnorm(5, 2)
      xi[1:4] ~ dCRP(1, size=4)
  })
  n <- 4
  data = list(y = rnorm(4))
  inits = list(xi = rep(1,4))
  model <- nimbleModel(code, data = data, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:4]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)

  
  ## various cases of clusters indep but not IID G0
  code <- nimbleCode({
      for(i in 1:3) {
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
      }
      thetaTilde[1] ~ dnorm(0, 1)
      thetaTilde[2] ~ dgamma(1, 1)
      thetaTilde[3] ~ dnorm(5, 1)
      xi[1:3] ~ dCRP(alpha, size = 3)
  })

  n <- 3
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:4]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)


  code <- nimbleCode({
      for(i in 1:3) {
          y[i] ~ dnorm(theta[i], 1)
          theta[i] <- thetaTilde[xi[i]]
      }
      thetaTilde[1] ~ dnorm(0, 1)
      thetaTilde[2] ~ dgamma(1, 1)
      thetaTilde[3] ~ dnorm(5, 1)
      xi[1:3] ~ dCRP(alpha, size = 3)
  })

  n <- 3
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:4]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  
  expect_identical(cn$numNodesPerCluster, 1L)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)


  ## case that is not allowed; doesn't like that thetas defined without using an indexing variable
  code <- nimbleCode({
      for(i in 1:3) {
          y[i] ~ dnorm(theta[i], 1)
          thetaTilde[i] ~ dnorm(0,1)
      }
      theta[1] <- thetaTilde[xi[1]]
      theta[2] <- exp(thetaTilde[xi[2]])
      theta[3] <- thetaTilde[xi[3]]
      xi[1:3] ~ dCRP(alpha, size = 3)
  })

  n <- 3
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Detected unusual indexing")

  ## another variation
  code <- nimbleCode({
      for(i in 1:4) {
          y[i] ~ dnorm(theta[i], 1)
          thetaTilde[i] ~ dnorm(0,1)
      }
      for(i in 1:2) 
          theta[i] <- thetaTilde[xi[i]]
      for(j in 3:4)
          theta[j] <- exp(thetaTilde[xi[j]])
      xi[1:4] ~ dCRP(alpha, size = 4)
  })

  n <- 4
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: differing number of clusters")

  ## cases of multiple obs per group regarding priors on cluster parameters

  ## clusters indep G0; multiple obs per group
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
              thetaTilde[i, j] ~ dnorm(i, 1)  # IID within cluster, indep across
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))


  ## clusters IID G0; indep but not identically distributed within cluster, 
  ## this uses conjugacy
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
              thetaTilde[i, j] ~ dnorm(j, 1)   # indep within cluster, IID across
          }}
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, each = J))


  ## clusters IID G0; indep but not identically distributed within cluster, 
  ## this case is of interest, discussed 2020-01-24 at Berkeley meeting

  code <- nimbleCode({
      for(i in 1:n) {
          y[i, 1] ~ dnorm(thetaTilde[xi[i], 1], 1)
          y[i, 2] ~ dgamma(thetaTilde[xi[i], 2], 1)
          thetaTilde[i, 1] ~ dnorm(0, 1)   # indep within cluster, IID across
          thetaTilde[i, 2] ~ dgamma(1, 1)   # indep within cluster, IID across
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 2
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = cbind(rnorm(n), rgamma(n,1,1)))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes[1:n])
  expect_identical(cn$clusterNodes[[2]], nodes[(n+1):(2*n)])
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))


  ## clusters IID G0, indep within cluster
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
          }
          thetaTilde[i, 1] ~ dnorm(0, 1)
          thetaTilde[i, 2] ~ dgamma(3, 1)
          thetaTilde[i, 3] ~ dnorm(5, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))

  
  ## clusters IID G0; dep within cluster
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
          }
          thetaTilde[i, 1] ~ dnorm(0,1)
          thetaTilde[i, 2] ~ dnorm(thetaTilde[i,1], 1)
          thetaTilde[i, 3] ~ dnorm(thetaTilde[i,2], 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, as.integer(J))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, J))

  ## mutilde and s2tilde; not conjugate
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i]])
              thetaTilde[i, j] ~ dnorm(0, 1)
          }
          s2tilde[i] ~ dinvgamma(1,1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = rgamma(n, 1, 1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, c(1L, 3L))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*(J+1)))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), c(1:n, rep(1:n, each = J)))


  ## mutilde and s2tilde; not conjugate, case 2
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i]], var = s2tilde[xi[i]])
          }
          thetaTilde[i] ~ dnorm(0, 1)
          s2tilde[i] ~ dinvgamma(1,1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n), s2tilde = rgamma(n, 1, 1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], nodes)
  expect_identical(cn$numNodesPerCluster, c(1L, 1L))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))

  ## mutilde and s2tilde; standard conjugacy                      
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i], j])
              thetaTilde[i, j] ~ dnorm(theta0, var = s2tilde[i, j]/kappa)
              s2tilde[i, j] ~ dinvgamma(1,1)
          }

      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = matrix(rgamma(J*n, 1, 1),n,J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, rep(as.integer(J), 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_invgamma_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = J), 2))



  ## mutilde and s2tilde; standard conjugacy but indexing is weird
  ## this has deps for s2tilde[1,1], but I think ok as changing s2tilde[1,1]
  ## doesn't impact likelihood so will amount to sample from prior
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i]+1, j])
              thetaTilde[i, j] ~ dnorm(theta0, var = s2tilde[i+1, j]/kappa)
          }
      }
      for(i in 1:(n+1)) {
          for(j in 1:J) {
              s2tilde[i, j] ~ dinvgamma(1,1)
          }
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = matrix(rgamma(J*(n+1), 1, 1),n+1,J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n+1, byrow = TRUE)[,-1]))
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, rep(as.integer(J), 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_invgamma_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))[-(1:J)]
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = J), 2))

  ## mutilde and s2tilde, IID across cluster, indep within, conjugate
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i], j])
              thetaTilde[i, j] ~ dnorm(j, var = s2tilde[i, j]/kappa)
              s2tilde[i, j] ~ dinvgamma(1,1)
          }

      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = matrix(rgamma(J*n, 1, 1),n,J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, rep(as.integer(J), 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_conjugate_dnorm_invgamma_dnorm")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = J), 2))

  ## mutilde and s2tilde, not IID across clusters
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i], j])
              thetaTilde[i, j] ~ dnorm(i, var = s2tilde[i, j]/kappa)
              s2tilde[i, j] ~ dinvgamma(1,1)
          }

      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = matrix(rgamma(J*n, 1, 1),n,J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], c(matrix(nodes, J, n, byrow = TRUE)))
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, rep(as.integer(J), 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(J))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*J*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(rep(1:n, each = J), 2))

  ## mutilde and s2tilde; not conj because don't have one parameter set per obs
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i]], var = s2tilde[xi[i]])
          }
          thetaTilde[i] ~ dnorm(theta0, var = s2tilde[i]/kappa)
          s2tilde[i] ~ dinvgamma(1,1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n), s2tilde = rgamma(n, 1, 1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(as.integer(1), 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*2))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), rep(1:n, 2))

  ## weird case of multiple thetaTilde but not s2tilde
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], var = s2tilde[xi[i]])
              thetaTilde[i,j] ~ dnorm(theta0, var = s2tilde[i]/kappa)
          }
          s2tilde[i] ~ dinvgamma(1,1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  J <- 3
  constants <- list(n = n, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = matrix(rnorm(J*n), n, J), s2tilde = rgamma(n, 1, 1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)

  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('s2tilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[2]], c(matrix(nodes, J, n, byrow = TRUE)))
  expect_identical(cn$numNodesPerCluster, c(1L, 3L))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, J)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers(c('thetaTilde', 's2tilde'))
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n*(J+1)))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), c(1:n, rep(1:n, each = J)))

  ## function in index - not allowed
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[n-i+1]], 1)
      }
      for(i in 1:n) {
          thetaTilde[i] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
                "findClusterNodes: Detected that a cluster parameter is indexed by a function")

  ## function in index - not allowed
  ## errors because can't figure out indexing of xi
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[exp(i)]], 1)
      }
      for(i in 1:n) {
          thetaTilde[i] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  expect_error(model <- nimbleModel(code, data = data, constants = constants, inits = inits),
               "dimensions specified are smaller than")


  ## function in index - not allowed
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[i]], s2tilde[xi[n-i+1]])
      }
      for(i in 1:n) {
          thetaTilde[i] ~ dnorm(0, 1)
      }
      for(i in 1:n)
          s2tilde[i] ~ dunif(0, 10)
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n), s2tilde = runif(n+1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_error(conf <- configureMCMC(model, print = FALSE),
               "findClusterNodes: Detected that a cluster parameter is indexed by a function")

  ## use of tilde var in two separate parameters; allowed, non-conj
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[i]], exp(thetaTilde[xi[i]]))
      }
      for(i in 1:n) {
          thetaTilde[i] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(1))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)

  ## use of tilde var in two separate parameters, partial intermediate; allowed, non-conj
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(theta[i], exp(thetaTilde[xi[i]]))
      }
      for(i in 1:n) {
          thetaTilde[i] ~ dnorm(0, 1)
          theta[i] <- thetaTilde[xi[i]]
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 1)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(2))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)


  ## case where inconsistent indexing would mess up cluster creation
  ## and nosample-empty-clusters
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[i]+1], exp(thetaTilde[xi[i]]))
      }
      for(i in 1:(n+1)) {
          thetaTilde[i] ~ dnorm(0, 1)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n+1))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Inconsistent indexing")

  ## two obs declarations so conj not detected
  code <- nimbleCode({
      for(i in 1:n) {
          y1[i] ~ dnorm(thetaTilde[xi[i]], 1)
          y2[i] ~ dnorm(thetaTilde[xi[i]], 1)
          thetaTilde[i] ~ dnorm(mu, sigma)
      }
      xi[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y1 = rnorm(n), y2 = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  cn <- nimble:::findClusterNodes(model, 'xi[1:5]')
  nodes <- model$expandNodeNames('thetaTilde')
  expect_identical(cn$clusterNodes[[1]], nodes)
  expect_identical(cn$numNodesPerCluster, rep(1L, 2))
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_silent(mcmc <- buildMCMC(conf))
  crpSampler <- mcmc$samplerFunctions[[crpIndex]]
  expect_equal(crpSampler$sampler, "CRP_nonconjugate")
  expect_identical(crpSampler$nObsPerClusID, 2)
  expect_identical(crpSampler$nIntermClusNodesPerClusID, as.integer(0))
  expect_identical(crpSampler$n, as.integer(n))

  paramSamplers <- conf$getSamplers('thetaTilde')
  expect_identical(sapply(paramSamplers, function(x) x$name), rep('CRP_cluster_wrapper', n))
  ids <- sapply(paramSamplers, function(x) x$control$clusterID)
  expect_identical(as.integer(ids), 1:n)


  ## thetaTildes indexed by multiple cluster variables - not allowed
  code <- nimbleCode({
      for(i in 1:n) {
          y[i] ~ dnorm(thetaTilde[xi[i]], 1)
          z[i] ~ dnorm(thetaTilde[eta[i]], 1)
          thetaTilde[i] ~ dnorm(0, 1)
      }                  
      xi[1:n] ~ dCRP(alpha, size = n)
      eta[1:n] ~ dCRP(alpha, size = n)
  })

  n <- 5
  constants <- list(n = n)
  data <- list(y = rnorm(n), z = rnorm(n))
  inits <- list(alpha = 1, xi = rep(1, n), eta = rep(1,n),
                thetaTilde = rnorm(n))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "sampler_CRP: Only the variables being clustered")

  ## too many xi values
  ## This error should be trapped more cleanly.
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
          }}
      for(i in 1:n) {
          for(j in 1:J) {
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      xi[1:(n+1)] ~ dCRP(alpha, size = n+1)
  })

  n <- 5
  n2 <- 4
  J <- 3
  constants <- list(n = n, n2 = n2, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n+1),
                thetaTilde = matrix(rnorm(J*n), n, J))
  model <- nimbleModel(code, data = data, constants = constants, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf), "replacement has length zero")

  ## too few xi values
  code <- nimbleCode({
      for(i in 1:n) {
          for(j in 1:J) {
              y[i, j] ~ dnorm(thetaTilde[xi[i], j], 1)
          }}
      for(i in 1:n) {
          for(j in 1:J) {
              thetaTilde[i, j] ~ dnorm(0, 1)
          }}
      xi[1:(n-1)] ~ dCRP(alpha, size = n-1)
  })

  n <- 5
  n2 <- 4
  J <- 3
  constants <- list(n = n, n2 = n2, J = J)
  data <- list(y = matrix(rnorm(n*J),n,J))
  inits <- list(alpha = 1, xi = rep(1, n-1),
                thetaTilde = matrix(rnorm(J*n), n, J))
  expect_error(model <- nimbleModel(code, data = data, constants = constants, inits = inits),
                                        "dimensions specified are smaller than model specification")


  ## Model 3 not yet available, but here are some initial cases

  ## Crossed clustering (model 3)

  code <- nimbleCode({
      for(i in 1:3) {
          for(j in 1:4) {
              y[i,j] ~ dnorm( thetaTilde[xi[i], eta[j]] , var = 1) 
              thetaTilde[i, j] ~ dnorm(0, 1)
          }
      }
      xi[1:3] ~ dCRP(1, size=3)
      eta[1:4] ~ dCRP(1, size=4)
  })
  inits <- list(xi = rep(1,3), eta = rep(1,4), 
                thetaTilde = matrix(0, nrow=3,  ncol=4))

  y <- matrix(5, nrow=3, ncol=4)
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf),
               "sampler_CRP: Detected use of multiple stochastic indexes of a variable")


  ## crossed clustering with truncation 
  code <- nimbleCode({
      for(i in 1:3) {
          for(j in 1:4) {
              y[i,j] ~ dnorm( thetaTilde[xi[i], eta[j]] , var = 1)
          }}
      for(i in 1:2) {
          for(j in 1:3) {
              thetaTilde[i, j] ~ dnorm(0, 1)
          }
      }
      xi[1:3] ~ dCRP(1, size=3)
      eta[1:4] ~ dCRP(1, size=4)
  })
  inits <- list(xi = rep(1,3), eta = rep(1,4), 
                thetaTilde = matrix(0, nrow=3,  ncol=4))

  y <- matrix(5, nrow=3, ncol=4)
  data <- list(y = y)
  model <- nimbleModel(code, data = data, inits = inits)
  expect_silent(conf <- configureMCMC(model, print = FALSE))
  expect_error(mcmc <- buildMCMC(conf),
               "sampler_CRP: Detected use of multiple stochastic indexes of a variable")


})



## simple tests of models

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dnorm(0,1)
    y[i] ~ dnorm(mu[xi[i]], sd = 1)
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rnorm(4))
data = list(y = rnorm(4))

testBUGSmodel(example = 'test1', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dpois(mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rgamma(4, 1, 1))
data = list(y = rpois(4, 4))

testBUGSmodel(example = 'test2', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dexp(mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rgamma(4, 1, 1))
data = list(y = rexp(4, 4))

testBUGSmodel(example = 'test3', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dgamma(4, mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rgamma(4, 1, 1))
data = list(y = rgamma(4, 4, 4))

testBUGSmodel(example = 'test4', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dbeta(1,1)
    y[i] ~ dbern(mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rbeta(4, 1, 1))
data = list(y = rbinom(4, size=1, prob=0.5))

testBUGSmodel(example = 'test5', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4){
    p[i,1:3] ~ ddirch(alpha=alpha0[1:3])
    y[i,1:3] ~ dmulti(prob=p[xi[i],1:3], size=3)
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
set.seed(1)
p0 <- matrix(0, ncol=3, nrow=4)
y0 <- matrix(0, ncol=3, nrow=4)
for(i in 1:4){
  p0[i,]=rdirch(1, c(1, 1, 1))
  y0[i,] = rmulti(1, prob=c(0.3,0.3,0.4), size=3)
}
inits = list(xi = 1:4, p=p0)
data = list(y = y0)
alpha0 = c(1,1,1)

testBUGSmodel(example = 'test6', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    s2[i] ~ dinvgamma(1, 1)
    mu[i] ~ dnorm(0,1)
    y[i] ~ dnorm(mu[xi[i]], var = s2[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=rnorm(4), s2=rinvgamma(4, 1,1))
data = list(y = rnorm(4))

testBUGSmodel(example = 'test7', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1, rate=1)
    y[i] ~ dinvgamma(shape=4, scale = mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = rep(1,4), mu=rgamma(4, 1, 1))
data = list(y = rinvgamma(4, 4, 4))

testBUGSmodel(example = 'test8', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dweib(shape=4, lambda = mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rgamma(4, 1, 1))
data = list(y = rweibull(4, 4, 4))

testBUGSmodel(example = 'test9', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dgamma(1,1)
    y[i] ~ dnorm(4, tau = mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rgamma(4, 1, 1))
data = list(y = rnorm(4, 4, 4))

testBUGSmodel(example = 'test10', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

model <- function() {
  for(i in 1:4) {
    mu[i] ~ dbeta(1,1)
    y[i] ~ dbinom(size=10, prob=mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rbeta(4, 1, 1))
data = list(y = rbinom(4, size=10, prob=0.5))

testBUGSmodel(example = 'test11', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)


model <- function() {
  for(i in 1:4) {
    mu[i] ~ dbeta(1,1)
    y[i] ~ dnegbin(size=10, prob=mu[xi[i]])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits =  list(xi = rep(1,4), mu=rbeta(4, 1, 1))
data = list(y = rnbinom(4, size=10, prob=0.5))

testBUGSmodel(example = 'test12', dir = "",
              model = model, data = data, inits = inits,
              useInits = TRUE)

model <- function() {
  for(i in 1:4){
    mu[i,1:4] ~ dmnorm(mu0[1:4], cov=Cov0[1:4, 1:4])
    y[i,1:4] ~ dmnorm(mu[xi[i],1:4], cov=Sigma0[1:4, 1:4])
  }
  xi[1:4] ~ dCRP(conc=1, size=4)
}
inits = list(xi = 1:4, mu=matrix(rnorm(16), 4, 4))
data = list(y = matrix(rnorm(16), 4, 4))
constants = list(mu0 = rep(0,4), Cov0 = diag(10, 4), Sigma0 = diag(1, 4))

testBUGSmodel(example = 'test13', dir = "",
              model = model, data = data, inits = c(inits, constants),
              useInits = TRUE)

## testing misspecification of dimension in a model

test_that("Testing of misspecification of dimension when using CRP", { 
  
  ## more labels than observations
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:10] ~ dCRP(conc=1, size=10)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,10), mu=rnorm(4)))
  conf <- configureMCMC(m)
  expect_error(buildMCMC(conf), "sampler_CRP: At least one variable has to be clustered")
  
  
  ## more observations than labels 
  code = nimbleCode({
    for(i in 1:10) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(conc=1, size=4)
  })
  expect_error(nimbleModel(code, data = list(y = rnorm(10)),
                           inits = list(xi = rep(1,4), mu=rnorm(10))),
               "dimensions specified are smaller")
  
  
  ## different obervations with same label
  code = nimbleCode({
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
    y[1] ~ dnorm(mu[xi[1]], 1)
    y[2] ~ dnorm(mu[xi[1]], 1)
    xi[1:2] ~ dCRP(conc=1, size=2)
  })
  m <- nimbleModel(code, data = list(y = rnorm(2)),
                   inits = list(xi = rep(1,2), mu=rnorm(2)))
  conf <- configureMCMC(m)
  expect_error(buildMCMC(conf), "sampler_CRP: Detected unusual indexing")
  
  
  ## same obervation with different label
  code = nimbleCode({
    mu[1] ~ dnorm(0,1)
    mu[2] ~ dnorm(0,1)
    y[1] ~ dnorm(mu[xi[1]], 1)
    y[1] ~ dnorm(mu[xi[2]], 1)
    xi[1:2] ~ dCRP(conc=1, size=2)
  })
  expect_error(nimbleModel(code, data = list(y = rnorm(2)),
                           inits = list(xi = rep(1,2), mu=rnorm(2))),
               "There are multiple definitions")
  
  ## less tilde variables than observations
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)  
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=1)  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50)))
  conf <- configureMCMC(m)
  expect_warning(buildMCMC(conf),
                 "sampler_CRP: The number of clusters based on the cluster parameters is less than the number of potential clusters")
  
  ## multiple tilde parameters
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)
      s2[i] ~ dinvgamma(1,1)
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=s2[xi[i]])  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50), s2=rinvgamma(50,1,1)))
  conf <- configureMCMC(m)
  expect_warning(buildMCMC(conf),
                 "sampler_CRP: The number of clusters based on the cluster parameters is less than the number of potential clusters")
  
  ## multiple tilde parameters, one is common for every observation
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)
      s2[i] ~ dinvgamma(1,1)
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=s2[xi[1]])  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50), s2=rinvgamma(1,1,1)))
  expect_error(conf <- configureMCMC(m), "findClusterNodes: found cluster membership parameters that use different indexing")
  
  ## more than one label used for each observation
  code = nimbleCode({
    for(i in 1:50){
      mu[i] ~ dnorm(0,1)  
    }
    for(i in 1:99){
      y[i] ~ dnorm(mu[xi[i]]+mu[xi[i+1]], var=1)  
    }
    y[100] ~ dnorm(mu[xi[100]], 1)
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = rnorm(100)),
                   inits = list(xi = rep(1,100), mu=rnorm(50)))
  expect_error(conf <- configureMCMC(m), "findClusterNodes: Detected that a cluster parameter is indexed by a function")
  
  ## test that a message is sent when truncation is hit
  code = nimbleCode({
    for(i in 1:3){
      mu[i] ~ dnorm(0,1)  
    }
    for(i in 1:100){
      y[i] ~ dnorm(mu[xi[i]], var=1)  
    }
    xi[1:100] ~ dCRP(conc=1, size=100)
  })
  m <- nimbleModel(code, data = list(y = c(rnorm(20, -5) , rnorm(20, 0), rnorm(20, 5),
                                           rnorm(20, 10), rnorm(20, 20))),
                   inits = list(xi = rep(1,100), mu=rnorm(3)))
  cm <- compileNimble(m)
  conf <- configureMCMC(m)
  expect_warning(mMCMC <- buildMCMC(conf))
  cmMCMC=compileNimble(mMCMC, project=m, resetFunctions=TRUE)
  set.seed(1)
  expect_output(cmMCMC$run(1), "CRP_sampler: This MCMC is for a parametric model")
  
})



## Test real BNP models:

test_that("Testing more BNP models based on CRP", { 
  ## Avandia meta-analysis
  codeBNP <- nimbleCode({
    for(i in 1:nStudies) {
      y[i] ~ dbin(size = nStudies, prob = q[i])
      x[i] ~ dbin(size = nStudies, prob = p[i])
      q[i] <- expit(theta + gamma[i])
      p[i] <- expit(gamma[i])
      gamma[i] ~ dnorm(mu[i], var = tau[i])
      mu[i] <- muTilde[xi[i]]
      tau[i] <- tauTilde[xi[i]]
    }
    for(i in 1:nStudies) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
      tauTilde[i] ~ dinvgamma(a0, b0)
    }
    xi[1:nStudies] ~ dCRP(conc, size = nStudies)
    conc ~ dgamma(1, 1)
    mu0 ~ dflat()
    sd0 ~ dunif(0, 100)
    a0 ~ dunif(0, 100)
    b0 ~ dunif(0, 100)
    theta ~ dflat()
  })
  
  Consts=list(nStudies=10)
  set.seed(1)
  Inits=list(gamma=rep(1,10),
             muTilde=rep(1,10),
             tauTilde=rep(1,10),
             xi=rep(1,10),
             conc =1,
             mu0 = 0,
             sd0 = 1,
             a0 = 1,
             b0 = 1,
             theta = 0)
  
  Data=list(y=rbinom(10, 10, 0.5), x=rbinom(10, 10, 0.5))
  
  model<-nimbleModel(codeBNP, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
  cmodel<-compileNimble(model)
  
  mConf <- configureMCMC(model)
  crpIndex <- which(sapply(mConf$getSamplers(), function(x) x[['name']]) == 'CRP')
  mMCMC <- buildMCMC(mConf)
  
  expect_equal(mConf$getSamplers()[[1]]$name, "CRP_concentration")
  expect_equal(class(mMCMC$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## Using testBUGSmodel
  model <- function() {
    for(i in 1:10) {
      y[i] ~ dbin(size = 10, prob = q[i])
      x[i] ~ dbin(size = 10, prob = p[i])
      q[i] <- expit(theta + gamma[i])
      p[i] <- expit(gamma[i])
      gamma[i] ~ dnorm(mu[i], var = tau[i])
      mu[i] <- muTilde[xi[i]]
      tau[i] <- tauTilde[xi[i]]
    }
    for(i in 1:10) {
      muTilde[i] ~ dnorm(mu0, sd = sd0)
      tauTilde[i] ~ dinvgamma(a0, b0)
    }
    xi[1:10] ~ dCRP(conc, size = 10)
    conc ~ dgamma(1, 1)
    mu0 ~ dflat()
    sd0 ~ dunif(0, 100)
    a0 ~ dunif(0, 100)
    b0 ~ dunif(0, 100)
    theta ~ dflat()
  }
  testBUGSmodel(example = 'test8', dir = "",
                model = model, data = Data, inits = Inits,
                useInits = TRUE)
  
  
  ## myeloma semiparametric AFT model
  ## data from 'myeloma' in the 'emplik' package
  ## library(emplik)
  ## data(myeloma)
  time <- c(1.25,1.25,2,2,2,3,5,5,6,6,6,6,7,7,7,9,11,11,11,11,11,13,14,15,16,16,17,17,18,19,19,24,25,26,32,35,37,41,41,51,52,54,58,66,67,88,89,92,4,4,7,7,8,12,11,12,13,16,19,19,28,41,53,57,77)
  vstatus <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0) # 0 = alive (i.e., censored)
  logBUN <- c(2.2175,1.9395,1.5185,1.7482,1.301,1.5441,2.2355,1.6812,1.3617,2.1139,1.1139,1.415,1.9777,1.0414,1.1761,1.7243,1.1139,1.2304,1.301,1.5682,1.0792,0.7782,1.3979,1.6021,1.3424,1.3222,1.2304,1.5911,1.4472,1.0792,1.2553,1.301,1,1.2304,1.3222,1.1139,1.6021,1,1.1461,1.5682,1,1.2553,1.2041,1.4472,1.3222,1.1761,1.3222,1.4314,1.9542,1.9243,1.1139,1.5315,1.0792,1.1461,1.6128,1.3979,1.6628,1.1461,1.3222,1.3222,1.2304,1.7559,1.1139,1.2553,1.0792)
  HGB <- c(9.4,12,9.8,11.3,5.1,6.7,10.1,6.5,9,10.2,9.7,10.4,9.5,5.1,11.4,8.2,14,12,13.2,7.5,9.6,5.5,14.6,10.6,9,8.8,10,11.2,7.5,14.4,7.5,14.6,12.4,11.2,10.6,7,11,10.2,5,7.7,10.1,9,12.1,6.6,12.8,10.6,14,11,10.2,10,12.4,10.2,9.9,11.6,14,8.8,4.9,13,13,10.8,7.3,12.8,12,12.5,14)
  
  n <- length(time)
  alive <- vstatus == 0
  cens_time <- rep(NA, n)
  cens_time[alive] <- time[alive]
  cens_time[!alive] <- Inf
  time[alive] <- NA
  
  logBUN <- (logBUN - mean(logBUN)) / sd(logBUN)
  HGB <- (HGB - mean(HGB)) / sd(HGB)
  
  ## accelerated failure time model per https://www4.stat.ncsu.edu/~ghosal/papers/PMR.pdf for Bayesian semiparametric AFT models
  codeAFT <- nimbleCode({
    for(i in 1:n) {
      x[i] ~ dweib(alpha, exp(lambda[i]))   # 'data' nodes
      is_cens[i] ~ dinterval(x[i], c[i])
      lambda[i] <-  inprod(Z[i, 1:p], delta[1:p]) + eta[i] 
      eta[i] <- etaTilde[xi[i]]
    }
    xi[1:n] ~ dCRP(conc, size = n)
    conc ~ dgamma(1, 1)
    for(i in 1:n){
      etaTilde[i] ~ dunif(b0, B0)
    }
    alpha ~ dunif(a0, A0)
    for(j in 1:p){
      delta[j] ~ dflat() 
    }
  })
  
  constants = list(b0 = -10, B0 = 10, a0 = 0.1, A0 = 10, p = 2, n = n, c
                   = cens_time, Z = cbind(logBUN, HGB))
  data = list(is_cens = as.numeric(alive), x = time)
  xInit <- rep(NA, n)
  xInit[alive] <- cens_time[alive] + 10
  inits = list(alpha = 1, delta = c(0, 0), conc = 1, 
               etaTilde = runif(n,constants$b0, constants$B0),
               xi = sample(1:3, n, replace = TRUE), x = xInit)
  
  model <- nimbleModel(codeAFT, constants = constants, data = data, inits = inits)
  conf = configureMCMC(model)
  mcmc = buildMCMC(conf)
  
  expect_equal(conf$getSamplers()[[1]]$name, "CRP_concentration")
  crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
  expect_identical(length(crpIndex), 1L)
  expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_nonconjugate")
  
  ## Using testBUGSmodel
  model <- function() {
    for(i in 1:n) {
      x[i] ~ dweib(alpha, 1+exp(lambda[i]))   # 'data' nodes
      is_cens[i] ~ dinterval(x[i], c[i])
      lambda[i] <-  inprod(Z[i, 1:p], delta[1:p]) + eta[i] 
      eta[i] <- etaTilde[xi[i]]
    }
    xi[1:n] ~ dCRP(conc, size = n)
    conc ~ dgamma(1, 1)
    for(i in 1:n){
      etaTilde[i] ~ dunif(b0, B0)
    }
    alpha ~ dunif(a0, A0)
    for(j in 1:p){
      delta[j] ~ dflat() 
    }
  }
  
  Data = list(is_cens = as.numeric(alive), x = time, 
              b0 = -10, B0 = 10, a0 = 0.1, A0 = 10, p = 2, n = n,
              c = cens_time, Z = cbind(logBUN, HGB))
  xInit <- rep(NA, n)
  xInit[alive] <- cens_time[alive] + 10
  Inits = list(alpha = 1, delta = c(0, 0), conc = 1, 
               etaTilde = runif(n,Data$b0, Data$B0),
               xi = sample(1:3, n, replace = TRUE), x = xInit)
  
  testBUGSmodel(example = 'test9', dir = "",
                model = model, data = Data, inits = Inits, 
                useInits = TRUE)
})


test_that("stick_breaking nimble function calculation and use is correct", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  ltruth <- log(truth)
  
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation"))
  
  expect_equal(stick_breaking(x, log=TRUE),
               ltruth,
               info = paste0("incorrect stick_breaking nimble function log calculation"))
  
  cSB <- compileNimble(stick_breaking)
  
  expect_equal(cSB(x, log=FALSE),
               truth,
               info = paste0("incorrect compiled stick_breaking nimble function calculation"))
  
  expect_equal(cSB(x, log=TRUE),
               ltruth,
               info = paste0("incorrect compiled stick_breaking nimble function log calculation"))
  
  x <- c(0.1, 0.4, -0.1, 0.3)
  expect_output(aux <- stick_breaking(x, log=FALSE), "values in 'z' have to be in", 
                info = "stick_breaking not warning of negative component")
  expect_equal(aux, rep(NaN, length(x)+1),
               info = "stick_breaking not correctly handling negative component")
  
  x <- c(0.1, 5, 0.4, 0.3)
  expect_output(aux <- stick_breaking(x, log=FALSE), "values in 'z' have to be in")
  expect_equal(aux, rep(NaN, length(x)+1),
               info = "stick_breaking incorrectly handling larger than 1 component")
  
  x <- c(0.1, 0.2, 0, 0.3, 0.8)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation with one 0 component"))
  
  x <- c(0.1, 0.2, 1, 0.3, 0.8)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  expect_equal(stick_breaking(x, log=FALSE),
               truth,
               info = paste0("incorrect stick_breaking nimble function calculation with one 1 component"))
})


test_that("Stick breaking model calculation is correct", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  
  SB_code <- nimbleCode({
    for(i in 1:5) z[i] ~ dbeta(1, 1)
    w[1:6] <- stick_breaking(z[1:5])
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  SB_model <- nimbleModel(SB_code, data=Inits)
  
  SB_model$z <- x
  SB_model$calculate()
  
  expect_equal(c(SB_model$w), truth,
               info = paste0("incorrect stick breaking weights in model"))
  
  c_SB_model <- compileNimble(SB_model)
  
  c_SB_model$z <- x
  c_SB_model$calculate()
  c_SB_model$w
  
  expect_equal(c(c_SB_model$w), truth,
               info = paste0("incorrect stick breaking weights in compiled model"))
  
})

## test simple models

model <- function() {
  for(j in 1:5) 
    z[j] ~ dbeta(1, 1)
  w[1:6] <- stick_breaking(z[1:5])
  for(i in 1:10){
    xi[i] ~ dcat(w[1:6])
  }
}

Inits <- list(z = rep(0.5,5))
Data <- list(xi = 1:10)

testBUGSmodel(example = 'test1', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)

model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    y[i] ~ dnorm( thetatilde[xi[i]], var=1)
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=1:10)
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test2', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)

model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    theta[i] <- thetatilde[xi[i]]
    y[i] ~ dnorm( theta[i], var=1)
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=1:10)
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test3', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)

model <- function(){
  for(i in 1:10){
    xi[i] ~ dcat(w[1:10])
    theta[i] <- thetatilde[xi[i]]
    y[i] ~ dnorm( theta[i], var=s2tilde[xi[i]])
  }
  for(i in 1:10){
    thetatilde[i] ~ dnorm(0, var=20)
    s2tilde[i] ~ dinvgamma(1, 1)
  }
  for(i in 1:9){
    z[i] ~ dbeta(1,1)
  }
  w[1:10] <- stick_breaking(z[1:9])
}

Inits=list(thetatilde=rep(0,10), z=rep(0.5, 9), xi=1:10, s2tilde=rep(1,10))
Data=list(y=rnorm(10))

testBUGSmodel(example = 'test4', dir = "",
              model = model, data = Data, inits = Inits,
              useInits = TRUE)



test_that("Testing conjugacy detection with bnp stick breaking models", { 
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dbeta(1, 1)
      }
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4)))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[6]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy")
  
  ## replace beta by uniform distribution: here no conjugacy is detected.
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dunif(0,1)
      }
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4)))
  conf <- configureMCMC(m)
  expect_failure(expect_match(conf$getSamplers()[[6]]$name,
                              "conjugate_dbeta_dcat",
                              info = "failed to detect categorical-beta conjugacy"))
  
  
  ## concentration parameter added in beta distribution
  code=nimbleCode(
    {
      for(i in 1:5){
        thetatilde[i] ~ dnorm(mean=0, var=10) 
      }
      for(i in 1:4){
        z[i] ~ dbeta(1, conc)
      }
      conc ~ dgamma(1,1)
      w[1:5] <- stick_breaking(z[1:4])
      
      for(i in 1:5){
        xi[i] ~ dcat(w[1:5])
        y[i] ~ dnorm(thetatilde[xi[i]], var=1)
      }
    }
  )
  m = nimbleModel(code, data = list(y = rnorm(5)),
                  inits = list(xi = rep(1,5), thetatilde=rep(0,5), z=rep(0.5,4), conc=1))
  conf <- configureMCMC(m)
  expect_match(conf$getSamplers()[[7]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy")
  
})


test_that("Testing BNP model using stick breaking representation", { 
  
  Code=nimbleCode(
    {
      for(i in 1:Trunc) {
        thetatilde[i] ~ dnorm(mean=0, var=40) 
        s2tilde[i] ~ dinvgamma(shape=1, scale=0.5) 
      }
      for(i in 1:(Trunc-1)) {
        z[i] ~ dbeta(1, 1)
      }
      w[1:Trunc] <- stick_breaking(z[1:(Trunc-1)])
      
      for(i in 1:N) {
        xi[i] ~ dcat(w[1:Trunc])
        theta[i] <- thetatilde[xi[i]]
        s2[i] <- s2tilde[xi[i]]
        y[i] ~ dnorm(theta[i], var=s2[i])
      }
    }
  )
  
  Consts <- list(N=50, Trunc=25)
  set.seed(1)
  Inits <- list(thetatilde = rnorm(Consts$Trunc, 0, sqrt(40)),
                s2tilde = rinvgamma(Consts$Trunc, shape=1, scale=0.5),
                z = rbeta(Consts$Trunc-1, 1, 1),
                xi = sample(1:10, size=Consts$N, replace=TRUE))
  Data = list(y = c(rnorm(Consts$N/2,5,sqrt(4)), rnorm(Consts$N/2,-5,sqrt(4))))
  
  model = nimbleModel(Code, data=Data, inits=Inits, constants=Consts,  calculate=TRUE)
  cmodel = compileNimble(model)
  
  modelConf = configureMCMC(model, thin=100)
  
  expect_match(modelConf$getSamplers()[[51]]$name, "conjugate_dbeta_dcat",
               info = "failed to detect categorical-beta conjugacy in BNP model")
  
  modelMCMC = buildMCMC(modelConf)
  CmodelMCMC = compileNimble(modelMCMC, project=model, resetFunctions=TRUE)
  
  CmodelMCMC$run(10000)
  
  ## results from the algorithm:
  samples = as.matrix(CmodelMCMC$mvSamples)
  s2Sam = samples[, 1:25]
  thetaSam = samples[, 26:50]
  zSam = samples[, 51:74]
  Tr = 25
  Wpost = t(apply(zSam, 1, function(x)c(x[1], x[2:(Tr-1)]*cumprod(1-x[1:(Tr-2)]), cumprod(1-x[1:(Tr-1)])[N=Tr-1])))
  
  ngrid = 302
  grid = seq(-10, 25,len=ngrid)
  
  ## posterior samples of the density
  nsave = 100
  predSB = matrix(0, ncol=ngrid, nrow=nsave)
  for(i in 1:nsave) {
    predSB[i, ] = sapply(1:ngrid, function(j)sum(Wpost[i, ]*dnorm(grid[j], thetaSam[i,],sqrt(s2Sam[i,]))))
  }
  
  ## hist(Data$y, freq=FALSE, xlim=c(min(grid), max(grid)))
  ## points(grid, apply(predSB, 2, mean), col="blue", type="l", lwd=2)
  
  f0 <- function(x) 0.5*dnorm(x,5,sqrt(4)) + 0.5*dnorm(x,-5,sqrt(4))
  fhat <- apply(predSB, 2, mean)
  f0grid <- sapply(grid, f0)
  
  L1dist <- mean(abs(f0grid - fhat))
  
  expect_equal(L1dist, 0.01, tol=0.01,
               info = "wrong estimation of density in DPM of normal distrbutions")
  
})



##-- test: random sampling from a compiled model adding one more level:
test_that("random sampling from model works fine", {
  set.seed(0)
  x <- rbeta(5, 1, 1)
  truth <- c(x[1], x[2:5]*cumprod(1-x[1:4]), prod(1-x[1:5]))
  
  SB_code2 <- nimbleCode({
    for(i in 1:5) 
      z[i] ~ dbeta(1, 1)
    w[1:6] <- stick_breaking(z[1:5])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  data <- list(xi = rep(1,10))
  SB_model2 <- nimbleModel(SB_code2, data=data, inits=Inits)
  
  c_SB_model2 <- compileNimble(SB_model2)
  
  c_SB_model2$z <- x
  c_SB_model2$calculate()
  
  expect_equal(c_SB_model2$w, truth,
               info = paste0("incorrect stick breaking weights in SB_model2"))
  
  #-- sampling via simulate:
  set.seed(0)
  simul_samp <- function(model) {
    model$simulate()
    return(model$w)
  }
  simul_samps <- t(replicate(10000, simul_samp(c_SB_model2)))
  
  trueE <- c(0.5^(1:5) )
  
  #-- checking the mean of the components of a vector that has a generalized dirichelt distribution
  #-- if z_i ~ beta then weights, w, defined by a SB representation have a generalized dirichlet distribution
  #-- and the expectation of w_j=0.5^j (in this case a=b=1).
  
  expect_equal(apply(simul_samps, 2, mean)[1:5], trueE, tol=0.01,
               info = paste0("incorrect weights (w) sampling  in SB_model2"))
  
  
  ## wrong specification of stick variables
  SB_code3 <- nimbleCode({
    for(i in 1:5) 
      z[i] ~ dgamma(10, 10)
    w[1:6] <- stick_breaking(z[1:5])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  ## wrong prior and starting values for stick variables
  set.seed(1)
  Inits <- list(z = rgamma(5, 10, 10))
  data <- list(xi = 1:10)
  expect_output(m <- nimbleModel(SB_code3, data=data, inits=Inits),
                "values in 'z' have to be in \\(0,1\\)")
  
  ## good starting values for stick variables
  set.seed(1)
  Inits <- list(z = rbeta(5, 1, 1))
  data <- list(xi = rep(1,10))
  SB_model3 <- nimbleModel(SB_code3, data=data, inits=Inits)
  expect_output(m <- SB_model3$simulate(), "values in 'z' have to be in \\(0,1\\)")
  
  ## wrong specification of length in stick variables, should be 5
  SB_code4 <- nimbleCode({
    for(i in 1:4) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:4])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  ## wrong length for stick variables
  set.seed(1)
  Inits <- list(z = rbeta(4, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(m <- nimbleModel(SB_code4, data=data, inits=Inits),
                 "number of items to replace")
  
  
  ## wrong specification of length in stick variables, should be 5
  ## no warning in nimbleModel function
  ## This is a shortcoming in NIMBLE processing of sizes in models.
  SB_code5 <- nimbleCode({
    for(i in 1:2) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:2])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  set.seed(1)
  Inits <- list(z = rbeta(2, 10, 10))
  data <- list(xi = rep(1,10))
  expect_failure(expect_error(SB_model5 <- nimbleModel(SB_code5, data=data, inits=Inits)))
  cSB_model5 <- compileNimble(SB_model5)
  expect_output(cSB_model5$calculate('w'), "Error in mapCopy")
  
  ## longer vector of stick variables
  SB_code6 <- nimbleCode({
    for(i in 1:10) 
      z[i] ~ dbeta(1,1)
    w[1:6] <- stick_breaking(z[1:10])
    for(i in 1:10){
      xi[i] ~ dcat(w[1:6])
    }
  })
  
  ## wrong length for stick variables, a warning is sent in the nimbleModel function
  set.seed(1)
  Inits <- list(z = rbeta(10, 10, 10))
  data <- list(xi = rep(1,10))
  expect_warning(SB_model6 <- nimbleModel(SB_code6, data=data, inits=Inits),
                 "number of items to replace")
  cSB_model6 <- compileNimble(SB_model6)
  expect_output(cSB_model6$calculate('w'), "Error in mapCopy")
})


## testing sampler assigment for conc parameter

test_that("Testing sampler assignment and misspecification of priors for conc parameter", { 
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dgamma(1, 1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  conf <- configureMCMC(m)
  expect_equal(conf$getSamplers('alpha')[[1]]$name, "CRP_concentration")
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dexp(1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  conf <- configureMCMC(m)
  expect_equal(conf$getSamplers('alpha')[[1]]$name, "RW")
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dunif(0,1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  conf <- configureMCMC(m)
  expect_equal(conf$getSamplers('alpha')[[1]]$name, "RW")
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(alpha, size=4)
    alpha ~ dnorm(-10,1) 
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4), alpha = 1))
  ## we do not warn of negative concentration values because there could be many such
  ## warnings in certain MCMC samplers for the concentration parameter
  expect_failure(expect_output(m$simulate(), "value of concentration parameter"))
  expect_output(out <- m$calculate(), "Warning: dynamic index out of bounds")
  ## think about better way to tell the user that the prior for alpha is wrong
  
  
  code = nimbleCode({
    for(i in 1:4) {
      mu[i] ~ dnorm(0,1)
      y[i] ~ dnorm(mu[xi[i]], sd = 1)
    }
    xi[1:4] ~ dCRP(0, size=4)
  })
  m = nimbleModel(code, data = list(y = rnorm(4)),
                  inits = list(xi = rep(1,4), mu = rnorm(4)))
  ## we do not warn of negative concentration values because there could be many such
  ## warnings in certain MCMC samplers for the concentration parameter
  expect_failure(expect_output(m$simulate(), "value of concentration parameter has to be larger than zero"))
  expect_output(out <- m$calculate(), "Warning: dynamic index out of bounds")
})


test_that("Testing dnorm_dnorm non-identity conjugacy setting, regression setting", { 

    ## Conjugacy detection and calculation of offset/coeff
    set.seed(1)
    code = nimbleCode({
        for(i in 1:n) {
            for(j in 1:J) {
                b1[i,j] ~ dnorm(beta, 0.25)
                y[i,j] ~ dnorm(b0 + b1[xi[i],j]*x[i], sd = 0.7)
            }
        }
        xi[1:n] ~ dCRP(conc=1, size=n)
        beta ~ dnorm(0,1)
        b0 ~ dnorm(0,1)
    })
    n <- 4
    J <- 2
    constants <- list(n = n, J = J)
    data <- list(y = matrix(rnorm(n*J),n,J), x = rnorm(n))
    m = nimbleModel(code, data = data, constants = constants,
                    inits = list(xi = c(4,3,2,1), b1 = matrix(rnorm(n*J),n,J), beta = rnorm(1), b0 = rnorm(1)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)
    crpIndex <- which(sapply(conf$getSamplers(), function(x) x[['name']]) == 'CRP')
    expect_equal(class(mcmc$samplerFunctions[[crpIndex]]$helperFunctions$contentsList[[1]])[1], "CRP_conjugate_dnorm_dnorm_nonidentity", info = 'dnorm_dnorm_nonidentity conjugacy not detected')
    mcmc$samplerFunctions[[1]]$helperFunctions[[1]]$calculate_offset_coeff(1,4)  # xi[1] = 4
    expect_identical(mcmc$samplerFunctions[[1]]$helperFunctions[[1]]$offset[1:J], rep(m$b0, J), info = 'calculation of offset in dnorm_dnorm_nonidentity incorrect')
    expect_equal(mcmc$samplerFunctions[[crpIndex]]$helperFunctions[[1]]$coeff[1:J], rep(m$x[1], J), tolerance = 1e-15, 
                 info = 'calculation of offset in dnorm_dnorm_nonidentity incorrect')
    mcmc$samplerFunctions[[1]]$helperFunctions[[1]]$calculate_offset_coeff(2,3)  # xi[2] = 3
    expect_identical(mcmc$samplerFunctions[[1]]$helperFunctions[[1]]$offset[1:J], rep(m$b0, J), info = 'calculation of offset in dnorm_dnorm_nonidentity incorrect')
    expect_equal(mcmc$samplerFunctions[[crpIndex]]$helperFunctions[[1]]$coeff[1:J], rep(m$x[2], J), tolerance = 1e-15,
                 info = 'calculation of offset in dnorm_dnorm_nonidentity incorrect')

    ## Correct predictive distribution
    tmp <- m$calculate()  ## in case we go back to having calculate_offset_coeff not recalculate after set to 0 and 1
    pYgivenT <- m$getLogProb('y[1, 1:2]')
    pT <- m$getLogProb('b1[4, 1:2]')  # 4 since xi[1]=4
    
    dataVar <- c(m$getParam('y[1, 1]', 'var'), m$getParam('y[1, 2]', 'var')) 
    priorVar <- c(m$getParam('b1[4, 1]', 'var'), m$getParam('b1[4, 2]', 'var'))
    priorMean <- c(m$getParam('b1[4, 1]', 'mean'), m$getParam('b1[4, 2]', 'mean'))
    postVar <- 1 / (m$x[1]^2 / dataVar + 1 / priorVar) # from conjugate sampler
    postMean <- postVar * (m$x[1]*(data$y[1, 1:2]-m$b0) / dataVar + priorMean / priorVar) # from conjugate sampler
    pTgivenY <- dnorm(m$b1[4, 1] , postMean[1], sqrt(postVar[1]), log = TRUE) +
                dnorm(m$b1[4, 2] , postMean[2], sqrt(postVar[2]), log = TRUE)  # from conjugate sampler

    mcmc$samplerFunctions[[1]]$helperFunctions$contentsList[[1]]$storeParams()
    mcmc$samplerFunctions[[1]]$helperFunctions$contentsList[[1]]$calculate_offset_coeff(1, 4)
    pY <- mcmc$samplerFunctions[[1]]$helperFunctions$contentsList[[1]]$calculate_prior_predictive(1)  
    expect_equal(pY, pT + pYgivenT - pTgivenY, info = "problem with predictive distribution for dnorm_dnorm_nonidentity")
    
    set.seed(1)
    mcmc$samplerFunctions[[1]]$helperFunctions$contentsList[[crpIndex]]$sample(1, 4)
    set.seed(1)
    smp <- rnorm(2, postMean, sqrt(postVar))
    expect_equal(smp, m$b1[4, 1:2], tolerance = 1e-15, info = "problem with predictive sample for dnorm_dnorm_nonidentity")

    ## Compare to identity conjugacy as special case.
    set.seed(1)
    n <- 100
    data <- list(y = matrix(rnorm(n*J),n,J), x = rnorm(n))
    constants <- list(n = n, J = J)
    inits <- list(xi = rep(1, n), b1 = matrix(4,n,J), beta = 1)
    
    set.seed(1)
    code = nimbleCode({
        for(i in 1:n) {
            for(j in 1:J) {
                b1[i,j] ~ dnorm(beta, 1)
                y[i,j] ~ dnorm(b0 + b1[xi[i],j]*x[i], sd = 1)
            }
        }
        xi[1:n] ~ dCRP(conc=1, size=n)
        beta ~ dnorm(0,1)
    })
    m = nimbleModel(code, data = data, constants = constants,
                    inits = c(inits, list(b0 = 0)))
    conf <- configureMCMC(m, monitors = c('b1','beta','xi'))
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    smp1 <- runMCMC(cmcmc, 1000, setSeed = 1)
    
    set.seed(1)
    code = nimbleCode({
        for(i in 1:n) {
            for(j in 1:J) {
                b1[i,j] ~ dnorm(beta, 1)
                y[i,j] ~ dnorm(b1[xi[i],j]*x[i], sd = 1)
            }
        }
        xi[1:n] ~ dCRP(conc=1, size=n)
        beta ~ dnorm(0,1)
    })
    m = nimbleModel(code, data = data, constants = constants,
                    inits = inits)
    conf <- configureMCMC(m, monitors = c('b1','beta','xi'))
    mcmc<-buildMCMC(conf)
    cm <- compileNimble(m)
    cmcmc <- compileNimble(mcmc, project = m)
    smp2 <- runMCMC(cmcmc, 1000, setSeed = 1)
    expect_identical(smp1, smp2, "sampling for identity and special case of non-identity not identical")
})    
   
 
test_that("Testing that cluster parameters are appropriately updated and mvSaved in good state", {
    ## Should always reject new clusters
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ T(dnorm(50, 1), -200, 200)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu = rep(0, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_identical(cm$mu[2], 0, info = 'mu[2] is changed')
    expect_identical(cm$mu[3], 0, info = 'mu[3] is changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())
    
    ## Should accept a cluster
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ T(dnorm(0, 1), -200, 200)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = c(0, rnorm(n-1, 50, 1)))
    inits <- list(xi = rep(1,n), mu = rep(50, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_equal(cm$mu[2], 0, tolerance = 3, info = 'mu[2] is not changed')
    expect_equal(cm$mu[3], 0, tolerance = 3, info = 'mu[3] is not changed')
    expect_identical(cm$mu[4], 50, info = 'mu[4] is changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())

    ## Should always reject new clusters
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(50, 1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu = rep(0, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_identical(cm$mu[2], 0, info = 'mu[2] is changed')
    expect_identical(cm$mu[3], 0, info = 'mu[3] is changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())
    
    ## Should accept a cluster
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0, 1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = c(0, rnorm(n-1, 50, 1)))
    inits <- list(xi = rep(1,n), mu = rep(50, n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(1)
    expect_equal(cm$mu[2], 0, tolerance = 3, info = 'mu[2] is not changed')
    expect_identical(cm$mu[3], 50, info = 'mu[3] is  changed')
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['logProb_mu']], cm$logProb_mu)
    expect_identical(sum(cm$logProb_mu), cm$calculate('mu'))
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(sum(cm$logProb_xi), cm$calculate('xi'))
    expect_equal(sum(cmcmc$mvSaved[['logProb_mu']]) + sum(cmcmc$mvSaved[['logProb_xi']]) + sum(cmcmc$mvSaved[['logProb_y']]), cm$calculate())

    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            muTilde[i] ~ T(dnorm(0, 1), -200, 200)
            sigmaTilde[i] ~ dinvgamma(a,b)
            mu[i] <- muTilde[xi[i]]
            y[i] ~ dnorm(mu[i], sd = sigmaTilde[xi[i]])
        }
        a ~ dgamma(1,1)
        b ~ dgamma(1,1)
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), muTilde = rnorm(n), sigmaTilde = rinvgamma(n, 1, 1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m, nodes = 'xi')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    cmcmc$run(10)
    expect_identical(cmcmc$mvSaved[['muTilde']], cm$muTilde)
    expect_identical(cmcmc$mvSaved[['logProb_muTilde']], cm$logProb_muTilde)
    expect_identical(cmcmc$mvSaved[['sigmaTilde']], cm$sigmaTilde)
    expect_identical(cmcmc$mvSaved[['logProb_sigmaTilde']], cm$logProb_sigmaTilde)
    expect_identical(cmcmc$mvSaved[['mu']], cm$mu)
    expect_identical(cmcmc$mvSaved[['xi']], cm$xi)
    expect_identical(cmcmc$mvSaved[['logProb_xi']], cm$logProb_xi)
    expect_identical(cmcmc$mvSaved[['y']], cm$y)
    expect_identical(cmcmc$mvSaved[['logProb_y']], cm$logProb_y)
    expect_identical(cmcmc$mvSaved[['a']], cm$a)
    expect_identical(cmcmc$mvSaved[['logProb_a']], cm$logProb_a)
    expect_identical(cmcmc$mvSaved[['b']], cm$b)
    expect_identical(cmcmc$mvSaved[['logProb_b']], cm$logProb_b)
  
})
  

test_that("Testing wrapper sampler that avoids sampling empty clusters", {
    set.seed(1)
    code = nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rnorm(n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[2]]$control$wrapped_type, 'conjugate_dnorm_dnorm_identity_dynamicDeps',
                     info = "cluster wrapper sampler not conjugate")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rgamma(n,1,1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[2]]$control$wrapped_type, 'RW',
                     info = "cluster wrapper sampler conjugate")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')


    ## Check that changes in cluster parameters persist even if no direct sampling of them
    
    set.seed(1)
    code = nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rnorm(n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    conf$removeSamplers('mu')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rgamma(n,1,1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    conf$removeSamplers('mu')
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (n+1):(2*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ']')
    focalClusterPresent <- apply(out[ , (n+1):(2*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    ## p=2 non-conjugate case
    set.seed(1)
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dnorm(0,1)
            sigma[i] ~ dinvgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], var = sigma[xi[i]])
        }
    })
    n <- 15
    data <- list(y = rnorm(n))
    inits <- list(xi = rep(1,n), mu = rnorm(n), sigma = rinvgamma(n, 1, 1))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[17]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ']')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (2*n+1):(3*n)])-1  # 2nd to last so have a few transitions
    focalClusterPresent <- apply(out[ , (2*n+1):(3*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    
    focalClusterName <- paste0('mu[', focalCluster, ']')
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')
    focalClusterName <- paste0('sigma[', focalCluster, ']')
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    ## ensure resetting cluster param values works for mv cluster parameter
    set.seed(1)
    code = nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i, 1:3] ~ dmnorm(z[1:3], pr[1:3,1:3])
            tmp[i, 1:3] <- exp(mu[xi[i],1:3])
            y[i,1:3] ~ dmnorm(tmp[i,1:3], pr[1:3,1:3])
        }
        
    })
    n <- 15
    data <- list(y = matrix(rnorm(n*3, 1, 1),n))
    inits <- list(xi = rep(1,n), mu = matrix(rnorm(n*3), n), pr = diag(3), z = rep(0,3))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    samplers <- conf$getSamplers()
    expect_identical(samplers[[2]]$name, 'CRP_cluster_wrapper',
                     info = "cluster wrapper sampler not set")
    expect_identical(samplers[[2]]$control$wrapped_type, 'RW_block',
                     info = "cluster wrapper sampler conjugate")
    cm <- compileNimble(m)
    mcmc <- buildMCMC(conf)
    cmcmc <- compileNimble(mcmc, project=m)
    out <- runMCMC(cmcmc, 500)
    expect_identical(1L, length(unique(out[ , paste0('mu[', n, ', 3]')])), info = "last cluster is sampled")
    focalCluster <- max(out[ , (3*n+1):(4*n)])-1  # 2nd to last so have a few transitions
    focalClusterName <- paste0('mu[', focalCluster, ', 3]')
    focalClusterPresent <- apply(out[ , (3*n+1):(4*n)], 1, function(x) focalCluster %in% x)
    focalClusterNew <- which(diff(focalClusterPresent) == 1)+1
    focalClusterAbsent <- which(!focalClusterPresent[-1])+1
    smp <- out[ , focalClusterName]
    expect_false(any(smp[focalClusterNew]- smp[focalClusterNew - 1] == 0),
                 info = 'no new cluster parameter value when new cluster opened')
    expect_true(all(smp[focalClusterAbsent] - smp[focalClusterAbsent - 1] == 0),
                info = 'new cluster parameter value despite cluster being closed')

    ## ensure that two different kinds of wrapped samplers compile and run
    code <- nimbleCode({
        xi[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu[i] ~ dgamma(1,1)
            y[i] ~ dnorm(mu[xi[i]], sd = 1)
        }
        xi2[1:n] ~ dCRP(conc=1, size=n)
        for(i in 1:n) {
            mu2[i] ~ dnorm(0,1)
            y2[i] ~ dnorm(mu2[xi2[i]], sd = 1)
        }
        
    })
    n <- 15
    data <- list(y = rnorm(n), y2 = rnorm(n))
    inits <- list(xi = rep(1,n), mu=rgamma(n,1,1), xi2 = rep(1,n), mu2 = rnorm(n))
    m <- nimbleModel(code, data=data, inits= inits, constants = list(n=n))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)
    cm <- compileNimble(m)
    ## It's hard to get testthat checking of logging to work right, so just
    ## running so that testthat catches if code errors out.
    cmcmc <- compileNimble(mcmc,project=m)
    out <- runMCMC(cmcmc, 10)
    
    
})


test_that("offset and coeff set up in conjugacy for BNP so that non-dependencies are screened out", {

    ## dnorm cases
    code <- nimbleCode({
        for(i in 1:2) {
            y[i] ~ dnorm(mu[xi[i]], 1)
            mu[i] ~ dnorm(0,1)
        }
        xi[1:2] ~ dCRP(1, 2)
    })
    m <- nimbleModel(code, data = list (y = rnorm(2)), 
                     inits = list(mu = rnorm(2), xi = rep(1,2)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)

    expect_identical(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]$N_dep_dnorm_identity, 2L)

    expect_identical(c('dep_dnorm_identity_offset', 'dep_dnorm_identity_coeff') %in%
                ls(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]), rep(TRUE, 2))

    ## dmnorm cases
    code <- nimbleCode({
        for(i in 1:2) {
            y[i, 1:3] ~ dmnorm(mu[xi[i], 1:3], pr[1:3,1:3])
            mu[i, 1:3] ~ dmnorm(z[1:3], pr[1:3,1:3])
        }
        xi[1:2] ~ dCRP(1, 2)
    })
    m <- nimbleModel(code, data = list (y = matrix(rnorm(6), 2)), 
                     inits = list(mu = matrix(rnorm(6),2), xi = rep(1,2), pr = diag(3)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)

    expect_identical(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]$N_dep_dmnorm_identity, 2L)

    expect_identical(c('dep_dmnorm_identity_offset', 'dep_dmnorm_identity_coeff') %in%
                ls(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]), rep(TRUE, 2))

    code <- nimbleCode({
        for(i in 1:2) {
            mn[i, 1:3] <- A[1:3,1:3]%*%mu[xi[i], 1:3]
            y[i, 1:3] ~ dmnorm(mn[i, 1:3], pr[1:3,1:3])
            mu[i, 1:3] ~ dmnorm(z[1:3], pr[1:3,1:3])
        }
        xi[1:2] ~ dCRP(1, 2)
    })
    m <- nimbleModel(code, data = list (y = matrix(rnorm(6), 2)), 
                     inits = list(A = diag(3), mu = matrix(rnorm(6),2), xi = rep(1,2), pr = diag(3)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)

    expect_identical(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]$N_dep_dmnorm_multiplicative, 2L)

    expect_identical(c('dep_dmnorm_multiplicative_offset', 'dep_dmnorm_multiplicative_coeff') %in%
                     ls(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]), rep(TRUE, 2))
    
    ## dwish cases
    code <- nimbleCode({
        for(i in 1:2) {
            y[i, 1:3] ~ dmnorm(mu[1:3], pr[xi[i], 1:3,1:3])
            pr[i, 1:3,1:3] ~ dwish(R[1:3,1:3], 8)
        }
        xi[1:2] ~ dCRP(1, 2)
    })
    pr <- array(0, c(2, 3, 3)); pr[1,,] <- pr[2,,] <- diag(3)
    m <- nimbleModel(code, data = list(y = matrix(rnorm(6),2)), 
                     inits = list(xi = rep(1,2), pr = pr, R = diag(3)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)

    expect_identical(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]$N_dep_dmnorm_identity, 2L)

    expect_identical(c('dep_dmnorm_identity_offset', 'dep_dmnorm_identity_coeff') %in%
                ls(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]), c(FALSE, TRUE))

    code <- nimbleCode({
        for(i in 1:2) {
            y[i, 1:3] ~ dmnorm(mu[1:3], pr0[i, 1:3, 1:3])
            pr0[i, 1:3,1:3] <- theta * pr[xi[i], 1:3,1:3]
            pr[i, 1:3,1:3] ~ dwish(R[1:3,1:3], 8)
        }
        xi[1:2] ~ dCRP(1, 2)
    })
    pr <- array(0, c(2, 3, 3)); pr[1,,] <- pr[2,,] <- diag(3)
    m <- nimbleModel(code, data = list(y = matrix(rnorm(6),2)), 
                     inits = list(xi = rep(1,2), pr = pr, R = diag(3)))
    conf <- configureMCMC(m)
    mcmc <- buildMCMC(conf)

    expect_identical(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]$N_dep_dmnorm_multiplicativeScalar, 2L)

    expect_identical(c('dep_dmnorm_multiplicativeScalar_offset', 'dep_dmnorm_multiplicativeScalar_coeff') %in%
                ls(mcmc$samplerFunctions[[2]]$regular_sampler[[1]]), c(FALSE, TRUE))
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
nimbleOptions(MCMCprogressBar = nimbleProgressBarSetting)

