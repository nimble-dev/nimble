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
  run = function(iiter = double(0), rangeii = double(0)){
    outVal <- mvSaved[tildeVar, iiter][rangeii]
    returnType(double(0))
    return(outVal)
  }
)


sampler_G2 <- nimbleFunction(
  # name = 'sampler_G'
  
  setup=function(model, mvSaved, useCompiled = TRUE){#, target, varNames
    # cheking the mvSaved object:
    m0 <- length(mvSaved$varNames) 
    if( m0 == 1 ){
      stop('you need at least one random variable depending on the random indexes') # 
    }
    
    # I need the names of the variables to access content in mvSaved.
    nodes <- model$getNodeNames() 
    VarNames <- model$getVarNames()
    
    detNodes <- model$getNodeNames(determOnly=TRUE) 
    stochNodes <- model$getNodeNames(stochOnly=TRUE) 
    dataNodes <- model$getNodeNames(dataOnly=TRUE) 
    distributions <- model$getDistribution(stochNodes) # finding  dCRP distr, if it exists, and the name of its node
    
    # getting variable and node xi: 
    # here I'm considering that xi is monitored and that there is only one dCRP distr in the model
    # -- there could be more than one dCRP distribution in the model. In this case we have more
    #    than one vector of xi, how do we identify the tilde variables for each?
    
    dCRPdistrindex <- distributions == 'dCRP'
    if(sum(dCRPdistrindex) == 1){
      dCRPDistr <- distributions[dCRPdistrindex]
      dCRPNode <- stochNodes[ dCRPdistrindex ] # xi nodes. Is the order in stochNodes and their distributions always the same?!?!?!!
      dCRPVar <- model$getVarNames(nodes=dCRPNode)
      #dCRPVarindex <- FALSE
      #i=1
      #while(dCRPVarindex == FALSE && i<=length(VarNames)){
      #  aux <- model$getDistribution(VarNames[i])
      #  aux2 <- aux[1] == 'dCRP'
      #  if(sum(aux2, na.rm=TRUE)==0){
      #    i=i+1
      #  }else{
      #    dCRPVar <- VarNames[i]
      #    dCRPVarindex <- TRUE
      #  }
      #}
    }else{
      if( sum(dCRPdistrindex) == 0 ){
        stop('there are no random indexes')
      }
      if( sum(dCRPdistrindex)>1 ){ # how to get tilde variables for each dCRP? Maybe no even a problem
        stop('only one dCRP distribution is allowed for now')
      }
    }
    
    
    
    # getting variable and node conc,  assuming there is only one conc parameter
    concNode <- FALSE
    i <- 1
    while(concNode == FALSE){ # finds the (first) node that depends on dCRPNode.
      aux <- model$getDependencies(nodes[i], self=FALSE)
      if(sum(aux == dCRPNode , na.rm=TRUE)==1){
        concNode <- nodes[i]
      }else{
        i <- i+1
      }
    } 
    concDistr <- model$getDistribution(concNode) # identifies if 'conc' is random or not
    if(is.na(concDistr)){
      concRnd <- FALSE
      concVar <- mvSaved$varNames[1] ## in this case, concVar won't be used, but must be a name in the mvSaved in order to compile
    }else{
      concRnd <- TRUE
      concVar <- concNode
    }
    # we could do this in te run code....
    AproxError <- 1e-10
    if(concRnd){ # defining truncation
      if(useCompiled){
        conc <- mvSaved$CobjectInterface[[concVar]]
      }else{
        conc <- mvSaved[[concVar]]
      }
      concHat <- mean(unlist(conc))
      Trunc <- ceiling(log(AproxError)/log(concHat/(concHat+1))+1) # gives an error of at most AproxError.
      algoConc <- 0
    }else{
      conc <- model$getParam(dCRPNode, 'conc') 
      Trunc <- ceiling(log(AproxError)/log(conc/(conc+1))+1)
      algoConc <- conc
    }
    
    N <- length(dataNodes) # sample size
    
    if(Trunc > N){
      print('for an approximation error smaller than 1e-10, Trunc > N.')
    }
    
  
    # getting tilde variables
    targetElements <- model$expandNodeNames(dCRPNode, returnScalarComponents=TRUE)
    #nInterm <- length(model$getDependencies(targetElements[1], determOnly = TRUE))
    tildeVarNames=c()
    itildeVar <- 1
    stochDepOnly <- FALSE
    
    Dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in 1:length(Dep)){ 
      Depi <- Dep[i]
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(Depi))
      expr <- parse(text = expr)[[1]]
      if(is.call(expr) && expr[[1]] == '[' && expr[[3]] == targetElements[1]){
        tildeVarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
        itildeVar <- itildeVar+1 
      }
    }
        
    
    #if(nInterm == 0){ stochDepOnly <- TRUE } # only stoch dependency on xi
    #if(nInterm >= 1){
    #  detDep <- model$getDependencies(targetElements[1], determOnly = TRUE)
    #  for(i in 1:length(detDep)){ # find the deterministic variables that depend xi
    #    detDepi <- detDep[i]
        #detDepExpr <- model$getValueExpr(detDepi)
    #    expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(detDepi))
    #    expr <- parse(text = expr)[[1]]
    #    if(is.call(expr) && expr[[1]] == '[' && expr[[3]] == targetElements[1]){
    #      tildeVarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
    #      itildeVar <- itildeVar+1 
    #    }
        #if(is.call(detDepExpr) && detDepExpr[[1]]=='['){ # deter. nodes that can depend on xi
        #  if(detDepExpr[[3]] == 'xi[1]'){
        #    tildeVarNames[itildeVar] <- VarNames[which(VarNames==detDepExpr[[2]])]
        #    itildeVar <- itildeVar+1 
        #  }
        #}
        #if(is.call(detDepExpr) && detDepExpr[[1]]!='[' && length(detDep)==1){ # another case of stochDepe only: where the deterministic node (only 1) is the reparam of variance 
        #  stochDepOnly <- TRUE # if length(detDep)>1 ????
        #}
    #  }
    #}
    #if(stochDepOnly){
    #  stochDep1 <- model$getDependencies(targetElements[1], self = FALSE) # y[1] and reparametrization
      #stochDep2 <- model$getDependencies(targetElements[1], self = FALSE, stochOnly=TRUE) # y[1]
    #  paramExprs1 <- nimble:::cc_getNodesInExpr(model$getValueExpr(stochDep1)) #e.g., 'mean = theta[xi[i]], var = sigma' 
      #paramExprs2 <- nimble:::cc_getNodesInExpr(model$getValueExpr(stochDep2)) #e.g., 'mean = theta[xi[i]], var = sigma' 
      
    #  for(i in 1:length(paramExprs1)){ # we get the tilde variance
    #    expr <- parse(text = paramExprs1[i])[[1]]
    #    if(is.call(expr) && expr[[1]]=='['){ # case where we have random sigma and its reparametrication; can it be >3?
    #      if(expr[[3]] == 'xi[1]'){
    #        tildeVarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
    #        itildeVar <- itildeVar+1 
    #      }
    #    }
    #  }
    #  for(i in 1:length(paramExprs2)){ # we get the tilde mean
    #    expr <- parse(text = paramExprs2[i])[[1]]
    #    if(is.call(expr) && expr[[1]]=='['){ # case where we have random sigm yand its reparametrication; can it be >3?
    #      if(expr[[3]] == 'xi[1]'){
    #        tildeVarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
    #        itildeVar <- itildeVar+1 
    #      }
    #    }
    #  }
    #}
    #tildeVarNames <- unique(tildeVarNames)
    
    p <- length(tildeVarNames)
    
    # storaging object:
    mvConf <- modelValuesConf(vars = 'G', type = 'double', size = list(G = c(Trunc, p+1)) )
    mv <- modelValues(mvConf, m = 1)
    getTildeVarList <- nimbleFunctionList(getTildeVarVirtual)
    #getTildeVarList <- nimble:::nimbleFunctionList(getTildeVarVirtual)
    for(j in 1:p){
      getTildeVarList[[j]] <- getTildeVar(mvSaved, tildeVarNames[j])
    }
  },
  
  run=function(){
    niter <- getsize(mvSaved)
    resize(mv, niter)
    for(iiter in 1:niter){
      if(concRnd==FALSE){
        conciter <- algoConc
      }else{
        conciter <- mvSaved[concVar, iiter][1]
      }
      #   
      #   #-- getting the unique values in the samples and their probabilities of beign sampled. Need for computing G later.
      probs <- numeric(N)
      uniqueValues <- matrix(0, ncol=p, nrow=N) # is this ok?
      xiiter <- mvSaved[dCRPVar, iiter]
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
      probs[index] <- conciter #probs <- probs/sum(probs)
      newvalueindex <- index
      
      #-- computing G:
      vaux <- rbeta(1, 1, conciter+N)
      v1prod <- 1
      Taux <- 1#Taux <- 0
      paramaux <- numeric(p)
      
      # first sampled value:
      index <- rcat(prob=probs[1:newvalueindex])
      if(index==newvalueindex){# sample from G_0
        model$simulate(tildeVarNames)
        for(j in 1:p){ 
          paramaux[j] <- values(model, tildeVarNames)[j]#values(model, tildeVarNames)[(j-1)*N+1]  #
        }
      }else{# sample one of the existing values
        for(j in 1:p){
          paramaux[j] <- uniqueValues[index, j] 
        }
      }
      for(j in 1:p){
        mv['G', iiter][Taux, j+1] <<- paramaux[j] # <<-
      }
      mv['G', iiter][Taux, 1] <<- vaux # <<-
      Taux <- Taux + 1
      
      # the rest of the values
      while(Taux <= Trunc-1){
        index <- rcat(prob=probs[1:newvalueindex])
        if(index==newvalueindex){# sample from G_0
          model$simulate(tildeVarNames)
          for(j in 1:p){ 
            paramaux[j] <- values(model, tildeVarNames)[j]#values(model, tildeVarNames)[(j-1)*N+1]  #
          }
        }else{# sample one of the existing values
          for(j in 1:p){
            paramaux[j] <- uniqueValues[index, j] 
          }
        }
        condaux <- mv['G', iiter][1:(Taux-1), 2] == paramaux[1]#uniqueValues[1:newvalueindex, 1] == paramaux[1] # check if we sample a new atom or an atom that ir in G already
        if(sum(condaux) >0){ # the atom already exists and we have to update the weights and not increase Trunc
          repindex=1
          while(condaux[repindex]==FALSE){
            repindex=repindex+1
          }
          v1prod <- v1prod*(1-vaux)
          vaux <-rbeta(1, 1, conciter+N)
          mv['G', iiter][repindex, 1] <<- mv['G', iiter][repindex, 1] + vaux*v1prod
        }else{ # agument the truncation and keep the same parameters
          for(j in 1:p){
            mv['G', iiter][Taux, j+1] <<- paramaux[j]
          }
          v1prod <- v1prod*(1-vaux)
          vaux <- rbeta(1, 1, conciter+N)
          mv['G', iiter][Taux, 1] <<- vaux*v1prod 
          Taux <- Taux+1
        }
      }
      # complete the vector of probabilities and atoms
      mv['G', iiter][Trunc, 1] <<- 1- sum(mv['G', iiter][1:(Trunc-1), 1])
      model$simulate(tildeVarNames)
      for(j in 1:p){ 
        mv['G', iiter][Trunc, j+1] <<- values(model, tildeVarNames)[j]#values(model, tildeVarNames)[(j-1)*N+1]  
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
    #conjugate <- FALSE # default
    
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents=TRUE)
    detDeps <- model$getDependencies(targetElements, determOnly = TRUE)
    
    # do this sampler only when conc ~ gamma(ac, bc) and only one dCRP distribution depends on conc
    #if(length(calcNodes) & length(detDeps)==0){ do this sampler!}else{M-H}
    #if(length(calcNodes)>2){ assign M-H sampler} # two dCRP distributions with same conc, or something like 'conc+conc1', conc1 random or not
    #if(length(detDeps)>0){assign M-H sampler} # conc param is deterministic
    
    aux <- targetElements==calcNodes
    xiNodes <- calcNodes[ aux == FALSE ]
    xiNodesi <- model$expandNodeNames(xiNodes, returnScalarComponents=TRUE)
    
    N <- length(xiNodesi) # number of observations
    conc <- model$getParam(xiNodes, 'conc')
    a <- model$getParam(target, 'shape')
    b <- model$getParam(target, 'rate')
    xi <- model[[xiNodes]]
    k <- length(unique(xi))
    
    ak <- a+k
    ak1 <- ak -1
    
    
  },
  
  
  run = function() {
    
    # -- generating augmented r.v. and computing the weight.
    x <- rbeta(1, conc+1, N)
    blog <- b-log(x)
    w <- ak1/(ak1 + N*blog)
    
    # -- updating the concentration parameter.
    if(runif(1)<=w){
      model[[target]] <<- rgamma(1, ak, blog)
    }else{
      model[[target]] <<- rgamma(1, ak1, blog)
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
        storeParams = function() {},
        calculate_prior_predictive = function(i = integer()) {
            returnType(double())
        },
        sample = function(i = integer()) {}
    )
)

CRP_nonconjugate <- nimbleFunction(
    contains = CRP_helper,
    setup = function(model, tildeVarNames, marginalizedNodes, dataNodes) {
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
            
## Claudia to create versions of this for the various conjugacies using code she has already
## written in the conjugacy-specific full samplers
CRP_conjugate_dnorm_dnorm <- nimbleFunction(
    contains = CRP_helper,
    setup = function(model, tildeVarNames, marginalizedNodes, dataNodes) {
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
            postVar <- 1/(1/dataVar + 1/priorVar)
            postMean <- postVar*(y/dataVar + priorMean/priorVar)
            model[[tildeVarNames]][i] <<- rnorm(1, postMean, sqrt(postVar)) # model[[marginalizedNodes[i]]] <<- rnorm(1, postMean, sqrt(postVar)) can not be accesed  when compiling the MCMC configuration
        }
    ))
        

CRP_conjugate_dgamma_dpois <- nimbleFunction(
  contains = CRP_helper,
  setup = function(model, tildeVarNames, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
    priora <- nimNumeric(1)
    priorb <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priora <<- model$getParam(marginalizedParam, 'shape') 
      priorb <<- model$getParam(marginalizedParam, 'rate') 
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      return(priora*log(priorb) - (priora+y)*log(priorb+1) + lgamma(priora+y) - lgamma(priora) - lfactorial(y))
    },
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])[1]
      model[[tildeVarNames]][i] <<- rgamma(1, shape=priora+y, rate=priorb+1)
    }
  ))


CRP_conjugate_dbeta_dbern <- nimbleFunction(
  contains = CRP_helper,
  setup = function(model, tildeVarNames, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
    priora <- nimNumeric(1)
    priorb <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priora <<- model$getParam(marginalizedParam, 'shape1') 
      priorb <<- model$getParam(marginalizedParam, 'shape2')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      return(lgamma(priora+y) + lgamma(priorb+1-y) - lgamma(priora) - lgamma(priorb) - log(priora+priorb))
    },
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])[1]
      model[[tildeVarNames]][i] <<- rbeta(1, shape1=priora+y, shape2=priorb+1-y)
    }
  ))


CRP_conjugate_dgamma_dexp <- nimbleFunction(
  contains = CRP_helper,
  setup = function(model, tildeVarNames, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
    priora <- nimNumeric(1)
    priorb <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priora <<- model$getParam(marginalizedParam, 'shape') 
      priorb <<- model$getParam(marginalizedParam, 'rate') 
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])[1]
      return(log(priora) + priora*log(priorb) - (priora+1)*log(priorb+y))
    },
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])[1]
      model[[tildeVarNames]][i] <<- rgamma(1, shape=priora+1, rate=priorb+y)
    }
  ))


CRP_conjugate_dgamma_dgamma <- nimbleFunction(
  contains = CRP_helper,
  setup = function(model, tildeVarNames, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
    priora <- nimNumeric(1)
    priorb <- nimNumeric(1)
  },
  methods = list(
    storeParams = function() {
      priora <<- model$getParam(marginalizedParam, 'shape') 
      priorb <<- model$getParam(marginalizedParam, 'rate')  
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      datashape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      return((datashape-1)*log(y) + priora*log(priorb) + lgamma(datashape+priora) -
               lgamma(datashape) - lgamma(priora) -(datashape+priora)*log(priorb+y))
    },
    sample = function(i = integer()) {
      datashape <- model$getParam(dataNodes[i], 'shape')
      y <- values(model, dataNodes[i])[1]
      model[[tildeVarNames]][i] <<- rgamma(1, shape=datashape+priora, rate=priorb+y)
    }
  ))


## Chris please check the following 2 comments
CRP_conjugate_ddirch_dmulti <- nimbleFunction(
  contains = CRP_helper,
  setup = function(model, tildeVarNames, marginalizedNodes, dataNodes) {
    ## this makes sure that we create objects to store persistent information used for all 'i'
    ## comment 1: check that this is the best way to get the dimension of the Dirichelt distribution's parameter
    d <- length(model[[marginalizedNodes[1]]])
    prioralpha <- nimNumeric(d)
  },
  methods = list(
    storeParams = function() {
      prioralpha <<- model$getParam(marginalizedParam, 'alpha')
    },
    calculate_prior_predictive = function(i = integer()) {
      returnType(double())
      y <- values(model, dataNodes[i])
      return(lfactorial(n) - sum(lfactorial(y)+lgamma(prioralpha))+
               lgamma(sum(prioralpha)) + sum(lgamma(prioralpha+y)) - lgamma(sum(prioralpha+y)))
    },
    sample = function(i = integer()) {
      y <- values(model, dataNodes[i])
      ## comment 2: check that this is the way to upadte the vector parameter
      model[[tildeVarNames]][i, ] <<- rdirch(alpha=prioralpha_y)
    }
  ))



# general dCRP sampler covering nonconjugate and conjugate cases
sampler_CRP <- nimbleFunction(
  name = 'sampler_CRP',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    ## note that even in inefficient case, we need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.
    
    ## browser()  # uncomment this to be able to step through code in debug mode ('n' for next line of code, 'f' to finish a loop or function, 'c' to continue to the next browser() call)

    ## Claudia - note that I've assumed essentially all the code below (up to the conjugacy checking I inserted) is the same for all of the different conjugate samplers as for the non-conjugate: yes, it is
      
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents=TRUE)
    n <- length(targetElements) # N2
    
    # first check that the sampler can be used: N=n, (n=N2)
    N <- length(model$getDependencies(targetElements, dataOnly = TRUE) )
    VarNames <- model$getVarNames()
    data <- model$getDependencies(targetElements[1], dataOnly = TRUE) 
    if(N != n){ stop('length of random indexes and observations has to be the same') }
    
    # finding tilde variables:
    tildeVarNames=c()
    itildeVar <- 1
    
    Dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in 1:length(Dep)){ 
      Depi <- Dep[i]
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(Depi))
      expr <- parse(text = expr)[[1]]
      ## see Chris' comments in email about reworking this a bit:
      #if(is.call(expr) && expr[[1]] == '[' && expr[[3]] == targetElements[1]){
      #  tildeVarNames[itildeVar] <- deparse(expr[[2]])#VarNames[which(VarNames==expr[[2]])] 
      #  itildeVar <- itildeVar+1 
      #}
      #if(is.call(expr) && expr[[1]] == '[' && all.vars(expr[[3]]) == model$getVarNames(nodes=target)){
      #    tildeVarNames[itildeVar] <- deparse(expr[[2]])
      #    itildeVar <- itildeVar+1 
      #}# all.vars(expr[[3]]) == model$getVarNames(nodes=target) does not work (logica(0)) when the models is reparametrized
      #for using the code proposed by Chris:
      cond1 <- all.vars(expr[[3]]) == model$getVarNames(nodes=target)
      if( length(cond) > 0  ){
        if( is.call(expr) && expr[[1]] == '[' ){
          tildeVarNames[itildeVar] <- deparse(expr[[2]])
          itildeVar <- itildeVar+1 
        }
      }
    }
    #tildeVarNames <- unique(tildeVarNames) # avoids repetition in reparametrized models
    

    
    N3 <- c()
    for(i in 1: length(tildeVarNames)){
      N3[i] <- length(model[[tildeVarNames[i]]])
    }
    
    if(sum(N3==N3[1]) != length(N3)){
      stop('tilde variables have different length')
    }
    for(i in 1:length(N3)){
      if(N3[i] < n){ stop('length of tilde node has to be at least equal to length of random indexes') }   
    }
    
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
          stop("Nimble cannot currently assign a sampler to a dCRP node unless each cluster indicator is associated with a single observation.")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
        if(length(detDeps) != nInterm) {
          type <- 'allCalcs'  # give up again; should only occur in strange situations
        } else {
          dataNodes[i] <- stochDeps[1]
          
          if(nInterm >= 1) {  # this should handle case of no intermediates - Chris
            intermNodes[i] <- detDeps[1]
            intermNodes2[i] <- detDeps[1]
            intermNodes3[i] <- detDeps[1]
          }
          if(nInterm >= 2)
            intermNodes2[i] <- detDeps[2]
          if(nInterm >= 3)
            intermNodes3[i] <- detDeps[3]
        }
      }
    }

      ## Claudia: please check that this correctly chooses nonconjugate if we have two tilde variables
      
      marginalizedNodes <- model$expandNodeNames(tildeVarNames[1])
      helperFunctions <- nimbleFunctionList(CRP_helper)
      # helperFunctions <- nimble:::nimbleFunctionList(CRP_helper)

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
      helperFunctions[[1]] <- eval(as.name(sampler))(model, tildeVarNames, marginalizedNodes, dataNodes)
      
      curLogProb <- numeric(n) # stores the los probabilities of sampling existing or not indicators
  },
  
  
  run = function() {
      ## Claudia, this assumes all code is the same for conjugate samplers and for nonconjugate except
      ## the three specific places 'helperFunctions[[1]]' is used. yes, it is
      conc <- model$getParam(target, 'conc')
      helperFunctions[[1]]$storeParams()
    #  -- udating xi:
    for(i in 1:n){ # updates one xi_i at the time , i=1,...,n
      xi_i <- model[[target]][i]
      xi <- model[[target]]
      cond <- sum(xi_i==xi) # if cond=1, xi_i is a singleton
      for(j in 1:n){ # calculate probability of sampling indexes 1,...,n   
        if(i==j){ # index i denotes a new indicator xi_i
          if(cond>1){ # a new parameter has to be created to calculate the prob
            newind <- 1
            mySum <- sum(xi == newind)
            while(mySum>0 & newind <= n) { # need to make sure don't go beyond length of vector
              newind <- newind+1
              mySum <- sum(xi == newind)
            }
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              if(nInterm >= 1) model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }else{ # we keep the old parameter as the "new" one
            newind <- xi_i
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              if(nInterm >= 1) model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes)  
          }
          curLogProb[j] <<- log(conc) + helperFunctions[[1]]$calculate_prior_predictive(i)
        }else{
          model[[target]][i] <<- model[[target]][j]
          if(type == 'indivCalcs') {
            if(nInterm >= 1) model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
          curLogProb[j] <<- model$getLogProb(dataNodes[i])
        }  
        model[[target]][i] <<- xi_i
      } # 
      
      index <- rcat(n=1, exp(curLogProb-max(curLogProb)))#
      if(index==i){# creates a new component: one that is not used
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

