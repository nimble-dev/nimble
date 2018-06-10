# BNP sampler: sampling labels in a DPM model when the random measure G is 
# integrated out. 
# Used when syntax xi[1:N] ~ dCRP(conc) is used in BUGS.

DP_measure = function( MCMCobject ) {
  model = MCMCobject$model  ## I think that the model object is accessible from somewhere in the MCMC object
  rsampler = nimble:::get_DP_measure_samples(model, MCMCobject$mvSamples)
  csampler = compileNimble(rsampler, project = model)
  csampler$run()
  samplesMeasure = csampler$samples
  
  namesVars <- rsampler$tildeVars
  p <- length(namesVars)
  truncG <- ncol(samplesMeasure) / (p+1)
  namesW <- sapply(1:truncG, function(i) paste( "weight[", i, "]", sep="" ))
  namesAtom <- unlist(sapply( 1:p, function(j) 
    sapply(1:truncG, function(i) paste( namesVars[j], "[", i, "]", sep="" )) ))
  
  colnames(samplesMeasure) <- c(namesW, namesAtom)
  return(samplesMeasure)
}



get_DP_measure_samples <- nimbleFunction(
  name = 'get_DP_measure_samples',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved){
    
    ## Check if the mvSaved is compiled or not.
    mvIsCompiled <- exists('dll', envir = mvSaved)
    if( mvIsCompiled ) {
      stop("get_DP_measure_samples: modelValues object has to be an uncompiled object.\n")
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
        stop('get_DP_measure_samples: One node with a dCRP distribution is required.\n')
      }
      stop('get_DP_measure_samples: Currently only models with one node with a dCRP distribution are allowed.\n')
    }
    if( sum(dcrpVar == mvSavedVars) == 0 ){
      stop( 'get_DP_measure_samples: The node having the dCRP distribution has to be monitored in the MCMC (and therefore stored in the modelValues object).\n')
    }
    
    
    ## Find the cluster variables, named tildeVars
    targetElements <- model$expandNodeNames(dcrpNode, returnScalarComponents = TRUE)
    tildeVars <- NULL
    itildeVar <- 1
    dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in seq_along(dep)) { 
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(dep[i])) 
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
    
    ## Check that cluster variables are monitored.
    counts <- sapply(tildeVars, function(x) x %in% mvSavedVars)  
    if( sum(counts) != length(tildeVars) ) {
      stop('get_DP_measure_samples: The node(s) representing the cluster variables has to be monitored in the MCMC (and therefore stored in the modelValues object).\n')  
    }
    
    if( is.null(tildeVars) ) { ## probably unnecessary as checked in CRP sampler, but best to be safe
      stop('get_DP_measure_samples: The model should have at least one cluster variable.\n')
    }
    
    ## Check that tilde variables are continuous and univariate variables (for avoiding long trials in simulating atoms for G in run code:
    for(i in seq_along(tildeVars)) {
      if( isDiscrete(model$getDistribution(tildeVars[i])[1]) ) { # the argument has to be the distribution of one node
        stop('get_DP_measure_samples: cluster variables should be continuous random variables.\n')
      }
      if( sum(model$isMultivariate(tildeVars)) > 0 ) {# getDimension(model$getDistribution(tildeVars[i])[1])
        stop( 'get_DP_measure_samples: only univariate cluster variables are allowed.\n' )
      }
    }
    
    ## Geting all parent nodes of cluster variables:
    ## Geting all parent nodes of cluster variables:
    parentNodesTildeVars <- NULL
    tildeVarsElements <- list()
    for(i in seq_along(tildeVars) ) {
      tildeVarsElements[[i]] <- model$expandNodeNames(tildeVars[i])
    }
    
    candidateParentNodes <- model$getNodeNames(includeData = FALSE)
    candidateParentNodes <- candidateParentNodes[!candidateParentNodes %in% unlist(tildeVarsElements)]
    for(i in seq_along(candidateParentNodes)) {
      aux <- model$getDependencies(candidateParentNodes[i], downstream = TRUE)
      for(j in seq_along(tildeVars)) {
        if( sum( aux == tildeVarsElements[[j]][1] ) )
          parentNodesTildeVars <- c(parentNodesTildeVars, candidateParentNodes[i])#stochNodes[i])
      }
    }
    ## Including cluster variables. Object to be used when simulating atoms in run code:
    ## Claudia note change, which ensures simulate() does its work in graphical order  
    parentNodesWithTildeVars <- model$expandNodeNames(c(parentNodesTildeVars, tildeVars), sort = TRUE)
    
    ## Getting  stochastic parent nodes of tilde variables used later for copying from mvSAved to model
    if( is.null(parentNodesTildeVars) ) {
      parentStochNodesTildeVars <- tildeVars
    } else {
      parentStochNodesTildeVars <- parentNodesTildeVars[model$isStoch(parentNodesTildeVars)]
      
      ## Checking that stochastic parent nodes of tilde variables are in mvSaved:
      counts <- sapply(parentStochNodesTildeVars, function(x) x %in% mvSavedVars)  
      if( sum(unlist(counts)) != length(parentStochNodesTildeVars) ) { # unlist in case of tilde variables with no parent nodes
        stop('get_DP_measure_samples: The stochastic parent nodes of the cluster variables have to be monitored in the MCMC (and therefore stored in the modelValues object).\n')  
      }
    }
    
    
    
    
    
    fixedConc <- TRUE # assume that conc parameter is fixed. This will change in the if statement if necessary
    
    ## Get all parents of xi, including deterministic cases such as theta[i] <- thetatilde[xi[i]]:
    parentNodes <- NULL
    nodesToCheck <- stochNodes[!stochNodes %in% c(dataNodes, dcrpNode)]  
    for(i in seq_along(nodesToCheck)){
      aux <- model$getDependencies(nodesToCheck[i], includeData = FALSE, stochOnly = TRUE)  
      if(sum(aux == dcrpNode)) 
        parentNodes <- c(parentNodes, aux[aux != dcrpNode])
    }
    
    if( length(parentNodes) ) { # concentration parameter is random
      ## Check that values stored in mvSaved are sufficient to calculate dcrp concentration.
      ## Do this by creating a model containing all NAs and then copying values from the mvSaved
      ## and then checking getParam gives a non-NA.
      fixedConc <- FALSE
      
      ## Which parent nodes are saved.
      savedParentNodes <- parentNodes[parentNodes %in% mvSavedVars]  
      
      ## Create model with NA values
      verbosity <- nimbleOptions('verbose')
      nimbleOptions(verbose = FALSE)
      modelWithNAs <- model$modelDef$newModel(check = FALSE, calculate = FALSE)
      nimbleOptions(verbose = verbosity)
      
      modelWithNAs[[dcrpNode]] <- model[[dcrpNode]] 
      ## copy savedParentNodes from mvSaved
      nimCopy(from = model, to = modelWithNAs, nodes = savedParentNodes) # Chris: should be mvSaved not model, we want to check that mvSaved has all we need to get conc, right? 
      if( length(savedParentNodes) == 0 ) { 
        stop( 'get_DP_measure_samples: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n') 
      } else {
        modelWithNAs$calculate(model$getDependencies(savedParentNodes)) 
        # Error in if (counts > 0) { : missing value where TRUE/FALSE needed
        dcrpParam <- modelWithNAs$getParam(dcrpNode, 'conc')
        if( is.na(dcrpParam) ) {
          stop('get_DP_measure_samples: Any variable involved in the definition of the concentration parameter must be monitored in the MCMC.\n')
        }
      }
      ## Determine stochastic and deterministic dependencies of parents of dCRP node
      allDepsOfSavedParentNodes <- model$getDependencies(savedParentNodes)
    } else { ## placeholder since allDepsOfSavedParentNodes must only have nodes in the mvSaved for correct compilation
      allDepsOfSavedParentNodes <- dcrpNode
      savedParentNodes <- dcrpNode # also used in run code
    }  
    
    N <- length(dataNodes)
    p <- length(tildeVars)
    nTilde <- length(values(model, tildeVars)) / p 
    # The error of approximation G wiht a is given by (conc / (conc +1))^{truncG-1}. 
    # we are going to define an error of aproximation and based on the posterior values of the conc parameter define the truncation level of G
    # the error is between errors that are considered very very small in the folowing papers
    # Ishwaran, H., & James, L. F. (2001). Gibbs sampling methods for stick-breaking priors. Journal of the American Statistical Association, 96(453), 161-173.
    # Ishwaran, H., & Zarepour, M. (2000). Markov chain Monte Carlo in approximate Dirichlet and beta two-parameter process hierarchical models. Biometrika, 87(2), 371-390.
    approxError <- 1e-10 
    
    ## Storage object to be sized in run code based on MCMC output (Claudia note change to comment)
    samples <- matrix(0, nrow = 1, ncol = 1)
    
    
  },
  
  run=function(){
    
    niter <- getsize(mvSaved) # number of iterations in the MCMC
    
    # defining the truncation level of the random measure's representation:
    if( fixedConc ) {
      dcrpAux <- model$getParam(dcrpNode, 'conc')
      concSamples <- rep(dcrpAux, niter)   ## Claudia, I would rename 'concSamples' as 'concSamples'
    } else {
      concSamples <- numeric(niter)
      for( iiter in 1:niter ) {
        nimCopy(from = mvSaved, to = model, nodes = savedParentNodes, row=iiter) 
        model$calculate(allDepsOfSavedParentNodes)
        concSamples[iiter] <- model$getParam(dcrpNode, 'conc')
      }
      dcrpAux <- mean(concSamples)
    }
    
    truncG <- log(approxError) / log(dcrpAux / (dcrpAux+1)) + 1
    truncG <- round(truncG)
    #approxError <- (dcrpAux / (dcrpAux +1))^(truncG-1)
    # I think is good to send message indicating what the truncation level is for an approximation error smaller than to 10^(-10)
    nimCat('get_DP_measure_samples: Approximating the random measure by a finite stick-breaking representation with and error smaller than 1e-10, leads to a truncation level of ', truncG, '.\n')
    
    ## Storage object: matrix with nrow = number of MCMC iterations, and ncol = (1 + p)*truncG, where
    ## truncG the truncation level of the random measure G (an integer given by the values of conc parameter)
    ## (p+1) denoted the number of parameters, each of length truncG, to be stored: p is the number of cluster components (length(tildeVars)) and 1 is for the weights.
    samples <<- matrix(0, nrow = niter, ncol = truncG*(p+1)) 
    
    
    for(iiter in 1:niter){
      
      ## getting the sampled unique values (tilde variables) and their probabilities of being sampled,
      ## need for computing density later.
      probs <- numeric(N)
      uniqueValues <- matrix(0, nrow = N, ncol = p)  
      xiiter <- mvSaved[dcrpVar, iiter]
      range <- min(xiiter):max(xiiter) 
      index <- 1
      for(i in seq_along(range)){   ## Chris changed to 'range' from 'rangei' - confusing because 'i' is used in multiple ways in the code here - 'iiter' and 'i'
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
      
      #-- computing G(.) = sum_{l=1}^{truncG} w_l delta_{atom_l} (.):
      vaux <- rbeta(1, 1, concSamples[iiter] + N)
      v1prod <- 1
      Taux <- 1
      paramAux <- numeric(p)
      
      #nimCopy(mvSaved, model, row = iiter) ## I think you need this - see comment below
      
      # first sampled values: w_1 and atom_1
      index <- rcat(prob = probs[1:newValueIndex])
      if(index == newValueIndex){   # sample from G_0
        
        ## Claudia if you simulate from the model then you are always using the current values of any parents of tildeVars, but don't you want to get the values from mvSaved for parents of the tildeVars? I think your test cases may always have the hyperparameters of the density of the tildevars be fixed, but I don't think this is always the case.
        ## Chris: now should be everithing ok with simulating from the model?
        ## Claudia - see comment about deterministic nodes in setup code 
        # we need to copy from mvSaved to model only once per iteration, right?
        nimCopy(mvSaved, model, parentStochNodesTildeVars, row = iiter) # copy posterior samples of parent nodes
        model$simulate(parentNodesWithTildeVars)
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
          model$simulate(parentNodesWithTildeVars)
          for(j in 1:p){ 
            paramAux[j] <- values(model, tildeVars)[(j-1)*nTilde + 1]
          }
        } else{  # sample one of the existing values
          for(j in 1:p){
            paramAux[j] <- uniqueValues[index, j] 
          }
        }
        condaux <- samples[iiter, truncG + 1:(Taux-1)] == paramAux[1]  # check if we sample a new atom or an atom that is in G already
        
        ## Claudia, if tildeVars are continuous, then if you sample from G_0 it will almost surely not equal to an existing atom, but if tildevars are  discrete variables, then if the sample from G_0 equals an existing value, is that considered a new atom or not? I guess not, in which case the code is fine, I think.
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
      
      # complete the vector of probabilities and atoms: w_Trunc and atom_Trunc
      samples[iiter, truncG] <<- 1 - sum(samples[iiter, 1:(truncG-1)])
      model$simulate(parentNodesWithTildeVars)
      for(j in 1:p){ 
        samples[iiter, (j+1)*truncG] <<- values(model, tildeVars)[(j-1)*nTilde + 1]
      }
    }
  },
  methods = list( reset = function () {} )
)




#-- Sampler for concentration parameter, conc, of the dCRP distribution.

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
    conc <- model[[target]] # model$getParam(xiNodes, 'conc')  
    xi <- model[[xiNode]]
    
    occupied <- numeric( N)
    for( i in 1:N )
      occupied[xi[i]] <- 1
    k <- sum(occupied)
    
    aux1 <- shapeParam + k
    aux2 <- aux1 - 1
    
    # -- generating augmented r.v. and computing the weight.
    x <- rbeta(1, conc+1, N)
    aux3 <- rateParam - log(x)
    w <- aux2/(aux2 + N*aux3)
    
    # -- updating the concentration parameter.
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
    sample = function(i = integer()) {}
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
    sample = function(i = integer()) {} ## nothing needed for non-conjugate
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
    sample = function(i = integer()) {
      dataVar <- model$getParam(dataNodes[i], 'var')
      y <- values(model, dataNodes[i])[1]
      postVar <- 1 / (1 / dataVar + 1 / priorVar)
      postMean <- postVar * (y / dataVar + priorMean / priorVar)
      ## Claudia, tildeVarNames and marginalizedNodes should be the same
      ## does this work and if so, I think we should remove 'tildeVarNames' in all of these conjugate cases
      model[[marginalizedVar]][i] <<- rnorm(1, postMean, sqrt(postVar)) # model[[marginalizedNodes[i]]] <<- rnorm(1, postMean, sqrt(postVar)) can not be accesed  when compiling the MCMC configuration
    }
  )
)

## Claudia: I changed parameter names to be more informative - is this ok and can you do it for all conjugacies?
## also added some spacing to the cod
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
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rgamma(1, shape = priorShape + y, rate = priorRate + 1)
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
    sample = function(i = integer()) {
      dataMean <- model$getParam(dataNodes[i], 'mean')
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rgamma(1, shape = priorShape + 0.5, rate = priorRate + (y-dataMean)^2/2)
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
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rbeta(1, shape1=priorShape1+y, shape2=priorShape2+1-y)
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
    sample = function(i = integer()) {
      dataSize <- model$getParam(dataNodes[i], 'size')
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rbeta(1, shape1=priorShape1+y, shape2=priorShape2+dataSize-y)
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
    sample = function(i = integer()) {
      dataSize <- model$getParam(dataNodes[i], 'size')
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rbeta(1, shape1=priorShape1+dataSize, shape2=priorShape2+y)
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
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rgamma(1, shape=priorShape+1, rate=priorRate+y)
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
    sample = function(i = integer()) {
      datashape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rgamma(1, shape=datashape+priorShape, rate=priorRate+y)
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
    sample = function(i = integer()) {
      dataShape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rgamma(1, shape=1+priorShape, rate=priorRate+y^dataShape)
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
    sample = function(i = integer()) {
      dataShape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      model[[marginalizedVar]][i] <<- rgamma(1, shape=dataShape+priorShape, rate=priorRate+1/y)
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
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])
      model[[marginalizedVar]][i, ] <<- rdirch(alpha=priorAlpha+y)
    }
  )
)



# general dCRP sampler covering nonconjugate and conjugate cases
sampler_CRP <- nimbleFunction(
  name = 'sampler_CRP',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    ## note that even in inefficient case, we need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.
    
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents=TRUE)
    targetVar <- model$getVarNames(nodes = target)  
    n <- length(targetElements) 
    
    # first check that the sampler can be used: we need one observation per random index
    nObs <- length(model$getDependencies(targetElements, stochOnly = TRUE, self = FALSE))
    if(n != nObs){ stop("sampler_CRP: The length of membership variable and observations has to be the same.\n") }
    
    ## finding 'tilde' variables (the parameters that are being clustered):
    tildeVars <- NULL
    itildeVar <- 1
    
    dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in seq_along(dep)) { 
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(dep[i])) 
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
      warning('sampler_CRP: The number of cluster parameters is less than the number of potential clusters. The MCMC is not strictly valid if ever it proposes more components than cluster parameters exist; NIMBLE will warn you if this occurs.\n')
    
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
        stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self=FALSE) 
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
    
    helperFunctions <- nimble:::nimbleFunctionList(CRP_helper)
    
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
                nimCat('CRP_sampler: This MCMC is not fully nonparametric. More components than cluster parameters exist are required.\n')
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
        ## Claudia: I don't think we need this next line - can you check?
        ## Answer: you are right, we don;t need it.
        #model[[target]][i] <<- xi[i]
      } # 
      
      index <- rcat(n=1, exp(curLogProb-max(curLogProb)))#
      if(index==i) {# creates a new component: one that is not used
          if(newind != xi[i]) {
              model[[target]][i] <<- newind
              helperFunctions[[1]]$sample(newind)
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

