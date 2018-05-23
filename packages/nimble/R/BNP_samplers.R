# BNP sampler: sampling labels in a DPM model when the random measure G is 
# integrated out. 
# Used when syntax xi[1:N] ~ dCRP(conc) is used in BUGS.

#-- we return a model Values object where each row has the sampled weights and atoms of measure G.
getTildeVarVirtual <- nimbleFunctionVirtual(
  run = function(iiter = double(0), rangeii = double(0))
    returnType(double(0))
)

getTildeVar <- nimbleFunction(
  contains = getTildeVarVirtual,
  setup = function(mvSaved, tildeVar){
  },
  run = function(iiter = double(0), rangeii = double(0)) {
    outVal <- mvSaved[tildeVar, iiter][rangeii]
    returnType(double(0))
    return(outVal)
  }
)


sampler_DP_density <- nimbleFunction(
  setup=function(model, mvSaved){
    
    ## Check if the mvSaved is compiled or not.
    mvIsCompiled <- exists('dll', envir = mvSaved)
    if( mvIsCompiled ) {
      stop("sampler_DP_density: modelValues object has to be an uncompiled object.\n")
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
        stop('sampler_DP_density: One node with a dCRP distribution is required.\n')
      }
      stop('sampler_DP_density: Currently only models with one node with a dCRP distribution are allowed.\n')
    }
    if( sum(dcrpVar == mvSavedVars) == 0 ){
      stop(paste('sampler_DP_density: The node having the dCRP distribution has to be monitored in the MCMC (and therefore stored in the modelValues object.\n'))
    }
    
    
    ## Find the 'tilde' variables (the parameters that are being clustered). 
    targetElements <- model$expandNodeNames(dcrpNode, returnScalarComponents = TRUE)
    tildeVars <- NULL
    itildeVar <- 1
    
    ## Claudia, shouldn't this code be the same as the block of code in CRP_sampler that we went back and forth on while refining it?
    ## also, here and in CRP_sampler, it's best to use 'seq_along(dep)' in case 'dep' is of zero length. Same for 1:length(expr) in CRP_sampler
    ## Answer: yes, here is the same code as in CRP_sampler using seq_along.
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
    
    ## Check that tilde variables are monitored.
    counts <- sapply(tildeVars, function(x) x %in% mvSavedVars)  ## Claudia, please check this looks ok. Answer: it does!
    if( sum(counts) != length(tildeVars) ) {
      stop(paste('sampler_DP_density: The nodes(s) representing the parameters that are being clustered must be monitored in the MCMC (and therefore stored in the modelValues object.\n'))  ## Claudia, please check this language in this message. 
    }
    if( is.null(tildeVars) ) {
      stop('sampler_DP_density: The model should have at least one tilde variable.\n')
    }
    
    fixedConc <- TRUE # assume that conc parameter is fixed. This will change in the if statement if necessary
    
    ## Determine stochastic and deterministic dependencies of parents of dCRP node  
    ## first get all dependencies of xi, including deterministic cases such as theta[i] <- thetatilde[xi[i]]: ## Claudia, aren't we getting deps of parents of xi not deps of xi? 
    allDep <- NULL
    for(i in seq_along(stochNodes)){
      aux <- model$getDependencies(stochNodes[i], includeData = FALSE, stochOnly = TRUE)  ## Claudia, ok, since dcrpNode would never be data?
      ## also, see below, can we simply add 'stochOnly = TRUE' above? Answer: yes, that is much shorter!
      if(sum(aux == dcrpNode)) { 
        ## aux <- model$getDependencies(stochNodes[i], includeData=FALSE)  ## ok to remove?
        allDep <- c(allDep, aux[aux != dcrpNode])
      }
    }
    
    ## Claudia, is this whole next block of code needed? In the code just above, can we do:
    ## aux <- model$getDependencies(stochNodes[i], includeData = FALSE, stochOnly = TRUE)
    ## actually, I'm confused about when we would even have 'theta' be in 'allDep'. Isn't allDep just going to be
    ## nodes that parameterize the dCRP concentration?
    
    # now eliminate deterministic dependencies ,e.g., theta[i] <- thetatilde[xi[i]], from allDep object: 
    #indexDetDep <- NULL
    #for(i in 1:length(tildeVars)) {  ## Claudia, why is this working with the tildeVars? Aren't we needing to work with dependencies of the parents of 'xi'?
    #  stochDepTildeVars <- model$getDependencies(tildeVars[i], determOnly=TRUE)
    #  if( sum(stochDepTildeVars[1] == allDep) >0 & is.na(sum(stochDepTildeVars[1] == allDep))==FALSE ) { # second condition checks that there are deterministic nodes
    #    for(j in 1:length(stochDepTildeVars)) {
    #      indexDetDep <- c(indexDetDep, which(stochDepTildeVars[j] == allDep))
    #    }
    #  }
    #}
    #if( is.null(indexDetDep) ) {
    #  parentNodes <- unique(allDep)  
    #} else {
    #  parentNodes <- unique(allDep[-indexDetDep])  
    #}
    parentNodes <- allDep
    
    if( length(parentNodes) ) { # concentration parameter not fixed
      ## Check that values stored in mvSaved are sufficient to calculate dcrp concentration.
      ## Do this by creating a model containing all NAs and then copying values from the mvSaved
      ## and then checking getParam gives a non-NA.
      fixedConc <- FALSE
      
      ## Which parent nodes are saved.
      ##savedParentNodes <- NULL
      ##for(i in seq_along(parentNodes)) {
      ##  savedParentNodes <- c(savedParentNodes, mvSavedVars[parentNodes[i]==mvSavedVars])
      ##}
      savedParentNodes <- parentNodes[parentNodes %in% mvSavedVars]  ## Claudia, ok to replace the above code with this? Answer: yes
      
      # model$simulate(savedParentNodes)  ## Claudia, why do we need this? Answer: we don't
      ## Create model with NA values
      verbosity <- nimbleOptions('verbose')
      nimbleOptions(verbose = FALSE)
      modelWithNAs <- model$modelDef$newModel(check = FALSE, calculate = FALSE)
      nimbleOptions(verbose = verbosity)
      
      #modelWithNAs[[dcrpNode]] <- model[[dcrpNode]] ## Claudia, why do we need this? Answer: we don't
      ## copy savedParentNodes from mvSaved
      nimCopy(from = model, to = modelWithNAs, nodes = savedParentNodes) ## Claudia, shouldn't this be savedParentNodes not parentNodes? Answer: Yes
      if( length(savedParentNodes) == 0 ) { 
        stop( paste('sampler_DP_density: Any variable involved in the concentration parameter must be monitored in the MCMC.\n') )
      } else {
        modelWithNAs$calculate(model$getDependencies(savedParentNodes)) 
        dcrpParam <- modelWithNAs$getParam(dcrpNode, 'conc')
        if( is.na(dcrpParam) ) {
          stop('sampler_DP_density: One or more variables related to concentration parameter are not being monitored in modelValues object.\n')
        }
      }
      allDepsOfSavedParentNodes <- model$getDependencies(savedParentNodes)
    } else allDepsOfSavedParentNodes <- dcrpNode  ## placeholder since allDepsOfSavedParentNodes must only have nodes in the mvSaved for correct compilation
    
    ## Claudia, should we check somewhere that there is at least one tildeVar in the model? Answer: Yes we should
    N <- length(dataNodes)
    p <- length(tildeVars)
    Ntilde <- length(values(model, tildeVars)) / p 
    approxError <- 1e-10 ## maximum allowable error in approximating unknown density with the truncation representation
    
    getTildeVarList <- nimble:::nimbleFunctionList(getTildeVarVirtual)
    for(j in 1:p){
      getTildeVarList[[j]] <- getTildeVar(mvSaved, tildeVars[j])
    }
    
    ## Storage object: matrix with nrow = number of MCMC iterations, and ncol = (1 + p)*Trunc, where
    ##    Trunc is an integer given by the values of conc parameter
    ##    size is a placeholder here as this will be set based on number of iterations in the mvSaved at run time
    samples <- matrix(0, nrow = 1, ncol = 1)
    
  },
  
  run=function(){
    
    niter <- getsize(mvSaved) # number of iterations in the MCMC
    
    #-- defining the Truncation level of the random measure's representation:
    # if 'conc' parameter id random, we get its value using steps c) to f)
    # if 'conc' is fixed, we get its value directly from the parameter of dCRP.
    # the relation between 'conc', trunc level, and error of approximation is: (conc / (conc +1))^{Trunc-1}=e 
    concSam <- numeric(niter)
    if( fixedConc == TRUE ) {
      concSam[1:niter] <- rep( model$getParam(dcrpNode, 'conc'), niter)
      dcrpAux <- model$getParam(dcrpNode, 'conc')
    } else {
      for( iiter in 1:niter ) {
        nimCopy(from = mvSaved, to = model, nodes = savedParentNodes, row=iiter) 
        # e): do calculate of model with NA
        model$calculate(allDepsOfSavedParentNodes)
        concSam[iiter] <- model$getParam(dcrpNode, 'conc')
      }
      dcrpAux <- mean(concSam)
    }
    
    Trunc <- log(AproxError)/log(dcrpAux/(dcrpAux+1)) + 1
    Trunc <- round(Trunc)
    
    setSize(samples, c(niter, Trunc*(p+1))) # first 1:Trunc columns are weights, then the atoms.
    
    for(iiter in 1:niter){
      #-- getting the sampled unique values (tilde variables) and their probabilities of beign sampled . Need for computing G later.
      probs <- numeric(N)
      uniqueValues <- matrix(0, ncol=p, nrow=N) # 
      xiiter <- mvSaved[dcrpVar, iiter]
      rangei <- min(xiiter):max(xiiter) # is this ok??
      index <- 1
      for(i in 1:length(rangei)){
        cond <- sum(xiiter==rangei[i])
        if(cond>0){
          probs[index] <- cond
          for(j in 1:p){
            uniqueValues[index, j] <- getTildeVarList[[j]]$run(iiter, rangei[i])
          }
          index <- index+1
        }
      }
      probs[index] <- concSam[iiter] #probs <- probs/sum(probs)
      newvalueindex <- index
      
      #-- computing G(.) = sum_{l=1}^{Trunc} w_l delta_{atom_l} (.):
      vaux <- rbeta(1, 1, concSam[iiter]+N)
      v1prod <- 1
      Taux <- 1
      paramaux <- numeric(p)
      
      # first sampled values: w_1 and atom_1
      index <- rcat(prob=probs[1:newvalueindex])
      if(index==newvalueindex){# sample from G_0
        model$simulate(tildeVars)
        for(j in 1:p){ 
          samples[ iiter, j*Trunc + Taux] <<- values(model, tildeVars)[(j-1)*Ntilde + 1]
        }
      }else{# sample one of the existing values
        for(j in 1:p){
          samples[ iiter, j*Trunc + Taux] <<- uniqueValues[index, j] 
        }
      }
      samples[ iiter, Taux] <<- vaux # <-
      Taux <- Taux + 1
      
      # the rest of the values: w_l and atom_l, l=1, ..Trunc-1
      while(Taux <= Trunc-1){
        index <- rcat(prob=probs[1:newvalueindex])
        if(index == newvalueindex){# sample from G_0
          model$simulate(tildeVars)
          for(j in 1:p){ 
            paramaux[j] <- values(model, tildeVars)[(j-1)*Ntilde + 1]
          }
        }else{# sample one of the existing values
          for(j in 1:p){
            paramaux[j] <- uniqueValues[index, j] 
          }
        }
        condaux <- samples[ iiter, Trunc + 1:(Taux-1)] == paramaux[1]#uniqueValues[1:newvalueindex, 1] == paramaux[1] # check if we sample a new atom or an atom that is in G already
        if(sum(condaux) >0){ # the atom already exists and we have to update the weights and not include a new value of the params 
          repindex=1
          while(condaux[repindex]==FALSE){
            repindex = repindex+1
          }
          v1prod <- v1prod*(1-vaux)
          vaux <-rbeta(1, 1, concSam[iiter]+N)
          samples[ iiter, repindex] <<- samples[iiter, repindex] + vaux*v1prod
        }else{ # agument the truncation and keep the same parameters
          for(j in 1:p){
            samples[ iiter, j*Trunc + Taux] <<- paramaux[j]
          }
          v1prod <- v1prod*(1-vaux)
          vaux <- rbeta(1, 1, concSam[iiter]+N)
          samples[iiter, Taux] <<- vaux*v1prod 
          Taux <- Taux+1
        }
      }
      
      # complete the vector of probabilities and atoms: w_Tunc and atom_Tunc
      samples[ iiter, Trunc] <<- 1 - sum(samples[ iiter, 1:(Trunc-1)])
      model$simulate(tildeVars)
      for(j in 1:p){ 
        samples[ iiter, (j+1)*Trunc] <<- values(model, tildeVars)[(j-1)*Ntilde + 1]#values(model, tildeVarNames)[(j-1)*N+1]  
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
    if(n != nObs){ stop("sampler_CRP: The length of random indexes and observations has to be the same.\n") }
    
    ## finding 'tilde' variables (the parameters that are being clustered):
    tildeVars <- NULL
    itildeVar <- 1
    
    dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in 1:length(dep)) { 
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(dep[i])) 
      for(j in 1:length(expr)) {
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

    nTilde <- sapply(tildeVars, function(x) length(model[[x]]))
    
    if(length(unique(nTilde)) != 1)
      stop('sampler_CRP: When multiple parameters are being clustered, the number of those parameters must all be the same.\n')

    min_nTilde <- min(nTilde) ## we need a scalar for use in run code
    if(min_nTilde < n)
      warning('sampler_CRP: The number of parameters to be clustered is less than the number of random indexes. The MCMC is not strictly valid if ever it proposes more components than exist; NIMBLE will warn you if this occurs.\n')
    
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
          stop("sampler_CRP: Nimble cannot currently assign a sampler to a dCRP node unless each cluster indicator is associated with a single observation.\n")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
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
                        conjugate_dgamma_dpois = 'CRP_conjugate_dgamma_dpois',
                        conjugate_dbeta_dbern  = 'CRP_conjugate_dbeta_dbern',
                        conjugate_dgamma_dexp = 'CRP_conjugate_dgamma_dexp',
                        conjugate_dgamma_dgamma = 'CRP_conjugate_dgamma_dgamma',
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
            while(mySum>0 & newind <= n) { # need to make sure don't go beyond length of vector
              newind <- newind+1
              mySum <- sum(xi == newind)
            }
            if(newind > min_nTilde) {
              nimCat('CRP_sampler: This MCMC is not fully nonparametric. More components than exist are required.\n')
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
        model[[target]][i] <<- newind
        helperFunctions[[1]]$sample(i)
      }else{
        model[[target]][i] <<- model[[target]][index]
      } 
    }
    model$calculate(calcNodes)
    
    
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {})
)

