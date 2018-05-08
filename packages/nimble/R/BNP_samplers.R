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
  setup=function(model, mvSaved, useCompiled = TRUE){#, target, varNames
    # cheking the mvSaved object:
    mvSavedVar <- mvSaved$varNames
    
    # We need the names of the variables to access content in mvSaved.
    nodes <- model$getNodeNames() 
    varNames <- model$getVarNames()
    
    stochNodes <- model$getNodeNames(stochOnly = TRUE)
    dataNodes <- model$getNodeNames(dataOnly = TRUE) 
    distributions <- model$getDistribution(stochNodes) # finding  dCRP distr, if it exists, and the name of its node
    
    #-- getting xi node and variable
    dcrpIndex <- which(distributions == 'dCRP')
    if(length(dcrpIndex) == 1) {
      dcrpNode <- stochNodes[dcrpIndex] 
      dcrpVar <- model$getVarNames(nodes = dcrpNode)
    } else {
      if(length(dcrpIndex) == 0 ){
        stop('sampler_DP_density: There are no random indexes. \n')
      }
      ## how to get tilde variables for each dCRP? Maybe not even a problem
      stop('sampler_DP_density: Currently only models with one node with a dCRP distribution are allowed. \n')
    }
    
    #-- checking that xi variable are monitored in mvSaved object:
    if( sum(dcrpVar == mvSavedVar) == 0 ){
      stop(paste('sampler_DP_density: labeling variable is not monitored in model values object. \n'))
    }
    
    
    #-- getting tilde variables
    targetElements <- model$expandNodeNames(dcrpNode, returnScalarComponents=TRUE)
    tildeVars <- NULL
    itildeVar <- 1
    
    dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in 1:length(dep)){ 
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(dep[i]))
      expr <- parse(text = expr)[[1]]
      foundTarget <- all.vars(expr[[3]]) == dcrpVar#targetVar# 
      if(length(foundTarget) > 0){
        if(is.call(expr) && expr[[1]] == '['){
          tildeVars[itildeVar] <- deparse(expr[[2]])
          itildeVar <- itildeVar+1 
        }
      }
    }
    
    #-- checking that tilde variables are monitored in mvSaved object:
    counts <- 0
    for(i in 1:length(tildeVars)) {
      if( sum(tildeVars[i] == mvSavedVar) == 1 ) {
        counts <- counts + 1
      }
    }
    if( counts != length(tildeVars) ) {
      stop(paste('sampler_DP_density: Some unique variable is not monitored in model values object. \n'))
    }
    

    #-- a) Stochastic and deterministic parents of xi
    # first get all dependencies of xi, including deterministic theta[i] <- thetatilde[xi[i]]:
    allDep <- c()
    for(i in 1:length(stochNodes)){
      aux <- model$getDependencies(stochNodes[i])
      if (sum(aux==dcrpNode) > 0) { 
        aux <- model$getDependencies(stochNodes[i], includeData=FALSE, )
        indexXi <- which(aux==dcrpNode)
        allDep <- c(allDep, aux[-indexXi])  
      }
    }
    if( length(allDep) == 0 ) { # conc parameter id fixed
      fixedConc <- TRUE
    } else { # conc parameter is random and need to find the parents of xi to check the mv object
      fixedConc <- FALSE
      # now eliminiate deterministic dependencies theta[i] <- thetatilde[xi[i]]: 
      indexDetDep <- c()
      for(i in 1:length(tildeVars)) {
        stochDepTildeVars <- model$getDependencies(tildeVars[i], determOnly=TRUE)
        if( sum(stochDepTildeVars[1] == allDep) >0 & is.na(sum(stochDepTildeVars[1] == allDep))==FALSE ) { # second condition checks that there are deterministic nodes
          for(j in 1:length(stochDepTildeVars)) {
            indexDetDep <- c(indexDetDep, which(stochDepTildeVars[j] == allDep))
          }
        }
      }
      if( is.null(indexDetDep) ) {
        parentNodes <- unique(allDep)  
      } else {
        parentNodes <- unique(allDep[-indexDetDep])  
      }
      
      #-- b) which parent nodes are in mvSaved:
      savedParentNodes <- c()
      for(i in 1:length(parentNodes)) {
        savedParentNodes <- c(savedParentNodes, mvSavedVar[parentNodes[i]==mvSavedVar])
      }
      model$simulate(savedParentNodes)
      # c) create model with NA values
      modelWithNAs <- model$modelDef$newModel(check = FALSE, calculate = FALSE)
      modelWithNAs[[dcrpNode]] <- model[[dcrpNode]] # if not included in NA model: Error in if (counts > 0) { : missing value where TRUE/FALSE needed
      nimCopy(from = model, to = modelWithNAs, nodes = parentNodes) # can nodes be a vector?
      # d) copy savedParentNodes from mvSaved: doing thi I get an error
      if( length(savedParentNodes) == 0 ) { # concentration parameter was not monitored
        stop( paste('sampler_DP_density: Some variable of conc parameter is not being monitored in model values object. \n') )
      } else {
        #for(i in 1:length(savedParentNodes)) {
          #value <- mvSaved[[savedParentNodes[i]]][[1]]
          #modelWithNAs[[savedParentNodes[i]]] <- value
          # e): do calculate of model with NA
          modelWithNAs$calculate(model$getDependencies(savedParentNodes)) #[i]
        #}
        # f)
        dcrpParam <- modelWithNAs$getParam(dcrpNode, 'conc')
        if( is.na(dcrpParam) ) {
          stop('sampler_DP_density: Some variable of conc parameter is not being monitored in model values object. \n')
        }
      }
    }
    
    allDepsOfSavedParentNodes <- model$getDependencies(savedParentNodes)

    N <- length(dataNodes)
    p <- length(tildeVars) 
    niter <- getsize(mvSaved) # number of iterations in the MCMC
    AproxError <- 1e-10 # error of approxiating g with the truncation represantation
    
    getTildeVarList <- nimble:::nimbleFunctionList(getTildeVarVirtual)
    for(j in 1:p){
      getTildeVarList[[j]] <- getTildeVar(mvSaved, tildeVars[j])
    }
    
    # storaging object: matrix with ncol = numer of iteration of MCMC and row = (1 + p)*Trunc, where
    #                   Trunc is an integer given by the values of conc parameter
    samples <- matrix(0, nrow=1, ncol=niter) # initial dimension of output matrix
  
    #mvConf <- modelValuesConf(vars = 'G', type = 'double', size = list(G = c(Trunc, p+1)) )
    #mv <- modelValues(mvConf, m = 1)
  
  },
  
  run=function(){
    
    #-- defining the Truncation level of the random measure's representation:
    # if 'conc' parameter id random, we get its value using steps c) to f)
    # if 'conc' is fixed, we get its value directly from the parameter of dCRP.
    # the relation between 'conc', trunc level, and error of approximation is: (conc / (conc +1))^{Trunc-1}=e 
    concSam <- numeric(niter)
    if( fixedConc ) {
      concSam[1:niter] <- nimNumeric( length = niter, value = model$getParam(dcrpNode, 'conc') )
      dcrpAux <- model$getParam(dcrpNode, 'conc')
    } else {
      for( iiter in 1:niter ) {
        for(i in 1:length(savedParentNodes)) {
          nimCopy(from = mvSaved, to = model, nodes = allDepsOfSavedParentNodes, row=iiter) # can nodes be a vector?
          #value <- mvSaved[[ savedParentNodes[i] ]][[iiter]]
          #model[[savedParentNodes[i]]] <- value
          # e): do calculate of model with NA
          model$calculate(model$getDependencies(savedParentNodes[i]))
        }
        concSam[iiter] <- model$getParam(dcrpNode, 'conc')
      }
      dcrpAux <- mean(concSam)
    }
    
    Trunc <- log(AproxError)/log(dcrpAux/(dcrpAux+1)) + 1
    Trunc <- round(Trunc)
    
    
    setSize(samples, c(Trunc*(p+1), niter)) # first 1:Trunc rows are weights, then the atoms.
    
    for(iiter in 1:niter){
      #-- getting the unique values in the samples and their probabilities of beign sampled. Need for computing G later.
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
      
      #-- computing G:
      vaux <- rbeta(1, 1, concSam[iiter]+N)
      v1prod <- 1
      Taux <- 1
      paramaux <- numeric(p)
      
      # first sampled value:
      index <- rcat(prob=probs[1:newvalueindex])
      if(index==newvalueindex){# sample from G_0
        model$simulate(tildeVars)
        for(j in 1:p){ 
          samples[j*Trunc + Taux, iiter] <- values(model, tildeVars)[j]
        }
      }else{# sample one of the existing values
        for(j in 1:p){
          samples[j*Trunc + Taux, iiter] <- uniqueValues[index, j] 
        }
      }
      samples[Taux, iiter] <<- vaux # <-
      Taux <- Taux + 1
      
      # the rest of the values
      while(Taux <= Trunc-1){
        index <- rcat(prob=probs[1:newvalueindex])
        if(index == newvalueindex){# sample from G_0
          model$simulate(tildeVars)
          for(j in 1:p){ 
            paramaux[j] <- values(model, tildeVars)[j]
          }
        }else{# sample one of the existing values
          for(j in 1:p){
            paramaux[j] <- uniqueValues[index, j] 
          }
        }
        condaux <- samples[Trunc + 1:(Taux-1), iiter] == paramaux[1]#uniqueValues[1:newvalueindex, 1] == paramaux[1] # check if we sample a new atom or an atom that is in G already
        if(sum(condaux) >0){ # the atom already exists and we have to update the weights and not include a new value of the params 
          repindex=1
          while(condaux[repindex]==FALSE){
            repindex = repindex+1
          }
          v1prod <- v1prod*(1-vaux)
          vaux <-rbeta(1, 1, concSam[iiter]+N)
          samples[repindex, iiter] <<- samples[repindex, iiter] + vaux*v1prod
        }else{ # agument the truncation and keep the same parameters
          for(j in 1:p){
            samples[j*Trunc + Taux, iiter] <<- paramaux[j]
          }
          v1prod <- v1prod*(1-vaux)
          vaux <- rbeta(1, 1, concSam[iiter]+N)
          samples[Taux, iiter] <<- vaux*v1prod 
          Taux <- Taux+1
        }
      }
      # complete the vector of probabilities and atoms
      samples[Trunc, iiter] <<- 1 - sum(samples[1:(Trunc-1), iiter])
      model$simulate(tildeVars)
      for(j in 1:p){ 
        samples[(j+1)*Trunc, iiter] <<- values(model, tildeVars)[j]#values(model, tildeVarNames)[(j-1)*N+1]  
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

