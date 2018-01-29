# BNP sampler: sampling labels in a DPM model when the random measure G is 
# integrated out. 
# Used when syntax xi[1:N] ~ dCRP(conc) is used in BUGS.

sampler_MarginalizedG_xi<-nimbleFunction(
  name = 'sampler_MarginalizedG_xi',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
    ## note that even in inefficient case, we need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.
    
    calcNodes <- model$getDependencies(target)
    targetElements <- model$expandNodeNames(target, returnScalarComponents=TRUE)
    n <- length(targetElements)
    
    dataNodes <- rep("", n)
    type <- 'indivCalcs'
    nInterm <- length(model$getDependencies(targetElements[1], determOnly = TRUE))
    intermNodes <- dataNodes
    intermNodes2 <- dataNodes
    intermNodes3 <- dataNodes
    if(nInterm > 3) type <- "allCalcs"  ## give up and do the inefficient approach
    for(i in seq_len(n)) {
      stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self=FALSE) 
      detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
      if(length(stochDeps) != 1) 
        stop("Nimble cannot currently assign a sampler to a dCRP node unless each cluster indicator is associated with a single observation.")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
      if(length(detDeps) != nInterm) {
        type <- 'allCalcs'  # give up again; should only occur in strange situations
      } else {
        dataNodes[i] <- stochDeps[1] 
        if(nInterm >= 1)
          intermNodes[i] <- detDeps[1]; intermNodes2[i] <- detDeps[1];
          intermNodes3[i] <- detDeps[1]
          if(nInterm >= 2)
            intermNodes2[i] <- detDeps[2]
          if(nInterm >= 3)
            intermNodes3[i] <- detDeps[3]
      }
    }
    
    curLogProb <- numeric(n) # stores the los probabilities of sampling existing or not indicators
    
  },
  
  
  run = function() {
    conc <- model$getParam(target, 'conc')
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
              model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }else{ # we keep the old parameter as the "new" one
            newind <- xi_i
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }
          curLogProb[j] <<- log(conc) + model$getLogProb(dataNodes[i]) #<<-
        }else{
          model[[target]][i] <<- model[[target]][j]
          if(type == 'indivCalcs') {
            model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
          curLogProb[j] <<- model$getLogProb(dataNodes[i]) #<<-
        }  
        model[[target]][i] <<- xi_i
      } # 
      
      index <- rcat(n=1, exp(curLogProb-max(curLogProb)))#
      if(index==i){# creates a new component: one that is not used
        model[[target]][i] <<- newind
      }else{
        model[[target]][i] <<- model[[target]][index]
      } 
    }
    model$calculate(calcNodes)
    
    
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {})
)

## setup code here determines if there is conjugacy in simple cases
## run code has not been changed, but will need to be added to to do the computation for the conjugate
## cases we decide to handle
## run code should be able to make use of 'conjugate', 'marginalizedParam' in 

sampler_MarginalizedG_general <- nimbleFunction(
  name = 'sampler_MarginalizedG_general',
  contains=sampler_BASE,
  
  setup=function(model, mvSaved, target, control){
      ## note that even in inefficient case, we need to do individual dataNodes[i] <- model$getDependencies(targetElements[i], stochOnly = TRUE) because we are not guaranteed that xi[i] is the cluster membership for y[i]; it could be xi[i] is associated with y[n-i+1], e.g.

      ## browser()  # uncomment this to be able to step through code in debug mode ('n' for next line of code, 'f' to finish a loop or function, 'c' to continue to the next browser() call)

      conjugate <- FALSE # default
      
      calcNodes <- model$getDependencies(target)
      targetElements <- model$expandNodeNames(target, returnScalarComponents=TRUE)
      n <- length(targetElements)
      nInterm <- length(model$getDependencies(targetElements[1], determOnly = TRUE))
      dataNodes <- rep("", n)
      type <- 'indivCalcs'
      
      intermNodes <- dataNodes
      intermNodes2 <- dataNodes
      intermNodes3 <- dataNodes
      if(nInterm > 3) type <- "allCalcs"  ## give up and do the inefficient approach
      for(i in seq_len(n)) {
          stochDeps <- model$getDependencies(targetElements[i], stochOnly = TRUE, self=FALSE) 
          detDeps <- model$getDependencies(targetElements[i], determOnly = TRUE)
          if(length(stochDeps) != 1) 
              stop("Nimble cannot currently assign a sampler to a dCRP node unless each cluster indicator is associated with a single observation.")  ## reason for this is that we do getLogProb(dataNodes[i]), which assumes a single stochastic dependent
          if(length(detDeps) != nInterm) {
              type <- 'allCalcs'  # give up again; should only occur in strange situations
          } else {
              dataNodes[i] <- stochDeps[1] 
              if(nInterm >= 1)
                  intermNodes[i] <- detDeps[1]; intermNodes2[i] <- detDeps[1];
              intermNodes3[i] <- detDeps[1]
              if(nInterm >= 2)
                  intermNodes2[i] <- detDeps[2]
              if(nInterm >= 3)
                  intermNodes3[i] <- detDeps[3]
          }
      }
 
      ## determination of conjugacy
    if(nInterm <= 1 && type == 'indivCalcs') {  # only detect conjugacy in simple settings for now
        conjugate <- TRUE  # at various steps below, may determine this not to be the case
        if(nInterm) {  # 1 intermediate; e.g. theta[i] <- thetatilde[xi[i]]
            detDep <- model$getDependencies(targetElements[1], determOnly = TRUE)  # e.g., 'theta[1]'
            detDepExpr <- model$getValueExpr(detDep)  # e.g., thetatilde[xi[1]]
            if(is.call(detDepExpr) && detDepExpr[[1]] == '[') {
                if(detDepExpr[[3]] != targetElements[1]) {
                    conjugate <- FALSE
                } else clusterNodes <- model$expandNodeNames(deparse(detDepExpr[[2]]))  # e.g., 'thetatilde[1]',...,
            }
        } else {  # 0 intermediates; e.g., y[i] ~ dnorm(theta[xi[i]], 1)
            stochDep <- model$getDependencies(targetElements[1], self = FALSE)  # e.g., y[1]
            paramExprs <- nimble:::cc_getNodesInExpr(model$getValueExpr(stochDep)) #e.g., 'mean = theta[xi[i]], var = sigma' 
            candidates <- lapply(paramExprs, function(expr) {
                expr <- parse(text = expr)[[1]]
                if(is.call(expr) && expr[[1]] == '[' && expr[[3]] == targetElements[1]) {
                    model$expandNodeNames(deparse(expr[[2]]))  # e.g., 'theta[1]', ...
                } else NULL
            })
            nonNull <- which(!sapply(candidates, is.null))
            if(length(nonNull) != 1) { # 'xi[1]' appears in multiple parameters of 'y[1]'
                conjugate <- FALSE
            } else clusterNodes <- candidates[[nonNull]]
        }
        if(conjugate) {  # not ruled out yet
            if(length(unique(model$getDeclID(clusterNodes))) != 1) { # cluster nodes are not all from same prior
                conjugate <-  FALSE
            } else {
                conjugacy <- model$checkConjugacy(clusterNodes[1], restrictLink = 'identity')
                ## 'restrictLink' is because we don't yet want to deal with cases like:
                ## y[i] ~ dnorm(5*theta[i], tau2)
                ## as the 5 complicates the marginal
                if(length(conjugacy)) {
                    conjugacy <- conjugacy[[1]]  # it's a one-element list of lists, so simplify it
                    conjugate <- TRUE
                    marginalizedParam <- conjugacy$target
                    ## print(conjugacy) # will show information you can use for conjugate sampler
                    ## I think you will simply want to use things like this in run code:
                    ## this would be the case for the simple normal-normal conjugacy
                    ## getParam(marginalizedParam, 'mean')
                    ## getParam(marginalizedParam, 'var') + getParam(dataNodes[i], 'var')
                    ## to construct the marginal you need
                } else conjugate <- FALSE
            }
        }
        cat("Conjugacy detected: ", conjugate)
    } 
      
    curLogProb <- numeric(n) # stores the los probabilities of sampling existing or not indicators
    
  },
  
  
  run = function() {
    conc <- model$getParam(target, 'conc')
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
              model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }else{ # we keep the old parameter as the "new" one
            newind <- xi_i
            model[[target]][i] <<- newind
            if(type == 'indivCalcs') {
              model$calculate(intermNodes[i])
              if(nInterm >= 2) model$calculate(intermNodes2[i])
              if(nInterm >= 3) model$calculate(intermNodes3[i])
              model$calculate(dataNodes[i])
            } else model$calculate(calcNodes) 
          }
          curLogProb[j] <<- log(conc) + model$getLogProb(dataNodes[i]) #<<-
        }else{
          model[[target]][i] <<- model[[target]][j]
          if(type == 'indivCalcs') {
            model$calculate(intermNodes[i])
            if(nInterm >= 2) model$calculate(intermNodes2[i])
            if(nInterm >= 3) model$calculate(intermNodes3[i])
            model$calculate(dataNodes[i])
          } else model$calculate(calcNodes) 
          curLogProb[j] <<- model$getLogProb(dataNodes[i]) #<<-
        }  
        model[[target]][i] <<- xi_i
      } # 
      
      index <- rcat(n=1, exp(curLogProb-max(curLogProb)))#
      if(index==i){# creates a new component: one that is not used
        model[[target]][i] <<- newind
      }else{
        model[[target]][i] <<- model[[target]][index]
      } 
    }
    model$calculate(calcNodes)
    
    
    copy(from = model, to = mvSaved, row = 1, nodes = calcNodes, logProb = TRUE)
  },
  methods = list( reset = function () {})
)


#-- Standarized output nimbleFunction
#-- user gives us: 
#   - model: the model
#   - mvSaved: modelValues object with the samples of xi, the tilde variables, and the conc parameter, if sampled 
#   - varNames: names of the variables in the mv object including the tilde variables and conc param (in that order!).
#   - rndconc: whether the concentration parameteris sampled or not (if not sampled is fixed)

#-- we return a model Values object (?) where each row is a matrix having the sampled weights and atom of measure G.

sampler_G <- nimbleFunction(
  # name = 'sampler_G'
  
  setup=function(model, mvSaved, target, varNames, rndconc){
    #calcNodes <- model$getDependencies(target) 
    targetElements <- model$expandNodeNames(target, returnScalarComponents=TRUE)
    
    # defining truncation for a fixed concentration parameter: 
    m0 <- length(varNames)
    if(rndconc==FALSE){ # fixed conc parameter
      conc <- model$getParam(target, 'conc') # works for deterministic and stochastic?
      AproxError <- 1e-10
      Trunc <- ceiling(log(AproxError)/log(conc/(conc+1))+1) # gives an error of at most AproxError.
    }else{
      conc <- mvSaved[[varNames[m0]]] #how to get the samples from this parameters? 
      concHat <- mean(unlist(conc))
      Trunc <- ceiling(log(AproxError)/log(concHat/(concHat+1))+1) # gives an error of at most AproxError.
    }
    
    niter <- length(mvSaved[[target]]) # number of iterations of the MCMC
    N <- length(targetElements) # sample size
    
    xi <- mvSaved[[target]]
    if(rndconc==FALSE){
      p <- length(mvSaved[[]])-1 #dimension of measure G when conc is not sampled
    }else{
      p <- length(mvSaved[[]])-2 #dimension of measure G when conc is  sampled
    }
    
    tildeVar1 <- rep(0, p)
    for(j in 1:p){
      tildeVar1[j] <- model$expandNodeNames(varNames[j])[1]
    }
    # storaging object:
    mvConf <- modelValuesConf(vars = 'G', type = 'double', size = list(G = c(Trunc, p+1)) )
    mv <- modelValues(mvConf, m = niter)
  },
  
  run=function(){
    for(iiter in 1:niter){
      
      if(rndconc==FALSE){
        conciter <- conc 
      }else{
        conciter <- mvSaved[[varNames[m0]]][[iiter]]
      }
      
      #-- getting the unique values in the samples and their probabilities of beign sampled. Need for computing G later.
      probs <- numeric(N)
      uniqueValues <- matrix(0, ncol=p, nrow=N) # is this ok?
      xiiter <- xi[[iiter]]
      rangei <- min(xiiter):max(xiiter) # is this ok??
      index <- 1
      for(i in 1:length(rangei)){
        cond <- sum(xiiter==rangei[i])
        if(cond>0){
          probs[index] <- cond
          for(j in 1:p){
            uniqueValues[index, j] <- mvSaved[[varNames[j]]][[iiter]][rangei[i]]
          }
          index <- index+1
        }
      }
      probs[index] <- conciter #probs <- probs/sum(probs)
      newvalueindex <- index
      
      #-- computing G: 
      vaux <- rbeta(1, 1, conciter+N)
      v1prod <- 1
      Taux <- 0
      paramaux <- numeric(p)
      while(Taux < Trunc-1){
        index <- rcat(prob=probs[1:newvalueindex])
        if(index==newvalueindex){# sample from G_0
          for(j in 1:p){ 
            model$simulate(tildeVar1[j])
            paramaux[j] <- values(model, varNames[j])[1]  
          }
        }else{# sample one of the existing values
          for(j in 1:p){
            paramaux[j] <- uniqueValues[index, j] 
          }
        }
        condaux <- uniqueValues[1:Taux, 1] == paramaux[1] # check if we sample a new atom or an atom that ir in G already
        if(sum(condaux) >0){ # the atom already exists and we have to update the weights and not increase Trunc
          repindex=1
          while(condaux[repindex]==FALSE){
            repindex=repindex+1
          }
          v1prod <- v1prod*(1-vaux)
          vaux <-rbeta(1, 1, conciter+N)
          mv['G', iiter][repindex, 1] <<- mv$G[[iiter]][repindex, 1] + vaux*v1prod
        }else{ # agument the truncation and keep the same parameters
          Taux <-Taux+1
          for(j in 1:p){
            mv$G[[iiter]][Taux, j+1] <<- paramaux[j]
          }
          if( Taux==1 ){
            mv['G', iiter][Taux, 1] <<-vaux
          }else{
            v1prod <- v1prod*(1-vaux)
            vaux <- rbeta(1, 1, conciter+N)
            mv['G', iiter][Taux, 1] <<-vaux*v1prod
          }
        }
      }
      # complete the vector of probabiities and atoms
      mv['G', iiter][Trunc, 1] <<- 1- sum(mv['G', iiter][1:(Trunc-1), 1])
      for(j in 1:p){ 
        model$simulate(varNames[j])
        mv['G', iiter][Trunc, j+1] <<- values(model, varNames[j])[1]  
      }
    }
  },
  
  methods = list( reset = function () {} )
)

#-- new veriosn of sampler_G where user only gives the model and the mv object
sampler_G2 <- nimbleFunction(
  # name = 'sampler_G'
  
  setup=function(model, mvSaved){#, target, varNames
    # cheking the mvSaved object:
    m0 <- length(mvSaved[[]]) 
    if( m0 == 1 ){
      stop('you need at least one random variable depending on the random indexes') # 
    }
    
    # I need the names of the variables to access content in mvSaved.
    nodes <- model$getNodeNames() # here are all the nodes in the model, in topological order.
    VarNames <- model$getVarNames()
    
    detNodes <- model$getNodeNames(determOnly=TRUE) 
    stochNodes <- model$getNodeNames(stochOnly=TRUE) 
    dataNodes <- model$getNodeNames(dataOnly=TRUE) 
    distributions <- model$getDistribution(stochNodes) # finding  dCRP distr, if it exists, and the name of its node
    
    # getting variable and node xi:
    dCRPdistrindex <- distributions == 'dCRP'
    if( sum(dCRPdistrindex) == 0 ){
      stop('there are no random indexes')
    }else{
      dCRPDistr <- distributions[dCRPdistrindex]
      dCRPNode <- stochNodes[ dCRPdistrindex ] # xi nodes. Is the order in stochNodes and their distributions always the same?!?!?!!
      dCRPVarindex <- FALSE
      i=1
      while(dCRPVarindex == FALSE){
        aux <- model$getDistribution(VarNames[i])
        aux2 <- aux[1] == 'dCRP'
        if(sum(aux2, na.rm=TRUE)==0){
          i=i+1
        }else{
          dCRPVar <- VarNames[i]
          dCRPVarindex <- TRUE
        }
      }
    }
    xi <- mvSaved[[dCRPVar]]
    
    # getting variable and node conc,  assuming there is only one conc parameter
    NodeDependXi <- c()
    for(i in 1:length(nodes)){ # might be inefficient
      aux <- model$getDependencies(nodes[i], self=FALSE)
      NodeDependXi[i] <- sum(aux == dCRPNode , na.rm=TRUE)
    }
    concNode <- nodes[NodeDependXi] 
    concDistr <- model$getDistribution(concNode)
    if(is.na(concDistr)){
      concRnd <- FALSE
    }else{
      concRnd <- TRUE
    }
    concVar <- concNode
    if(concRnd){ # defining truncation
      conc <- mvSaved[[concVar]]
      concHat <- mean(unlist(conc))
      Trunc <- ceiling(log(AproxError)/log(concHat/(concHat+1))+1) # gives an error of at most AproxError.
    }else{
      conc <- model$getParam(dCRPNode, 'conc') 
      AproxError <- 1e-10
      Trunc <- ceiling(log(AproxError)/log(conc/(conc+1))+1)
    }
    
    targetElements <- model$expandNodeNames(dCRPNode, returnScalarComponents=TRUE)
    niter <- length(xi) # number of iterations of the MCMC
    N <- length(targetElements) # sample size
    if(concRnd==FALSE){
      p <- length(mvSaved[[]])-1 #dimension of measure G when conc is not sampled
    }else{
      p <- length(mvSaved[[]])-2 #dimension of measure G when conc is  sampled
    }
    
    # getting tilde variables by finding random nodes with N components
    tildevarNames=c()
    i=1
    for(j in 1:length(VarNames)){
      aux <- model$getDistribution(VarNames[j]) # can 
      if(length(aux) == N){ # VarNames[j] has N components
        if(is.na(aux[1])==FALSE){ # VarNames[j] is a random variable 
          aux2 <- model$expandNodeNames(VarNames[j])
          if(aux2[1]!=dataNodes[1]){  # VarNames[j] is different from data
            tildevarNames[i]=VarNames[j]
            i=i+1
          }
        }
      }
    }
    # there could be other random nodes with N components
    if(length(tildevarNames)!=p){ # matching the names in tildevarNames with anme sin mvSaved
      cat('warning: there might be some ERRORs that can be safely solved: have to TEST!')
      tildevarNames2=c()
      i=1
      for(j in 1:length(tildevarNames)){
        if(length(mvSaved[[tildevarNames[j]]])==niter){
          tildevarNames2[i]=tildevarNames[j]
          i=i+1
        }
      }
      if(length(tildevarNames2)!=p){
        cat('warnings: do not know what is happening??')
      }
    }
    
    
    # storaging object:
    mvConf <- modelValuesConf(vars = 'G', type = 'double', size = list(G = c(Trunc, p+1)) )
    mv <- modelValues(mvConf, m = niter)
  },
  
  run=function(){
    for(iiter in 1:niter){
      
      if(concRnd==FALSE){
        conciter <- conc 
      }else{
        conciter <- conc[[iiter]]
      }
      
      #-- getting the unique values in the samples and their probabilities of beign sampled. Need for computing G later.
      probs <- numeric(N)
      uniqueValues <- matrix(0, ncol=p, nrow=N) # is this ok?
      xiiter <- xi[[iiter]]
      rangei <- min(xiiter):max(xiiter) # is this ok??
      index <- 1
      for(i in 1:length(rangei)){
        cond <- sum(xiiter==rangei[i])
        if(cond>0){
          probs[index] <- cond
          for(j in 1:p){
            uniqueValues[index, j] <- mvSaved[[tildevarNames[j]]][[iiter]][rangei[i]]
          }
          index <- index+1
        }
      }
      probs[index] <- conciter #probs <- probs/sum(probs)
      newvalueindex <- index
      
      #-- computing G: 
      vaux <- rbeta(1, 1, conciter+N)
      v1prod <- 1
      Taux <- 0
      paramaux <- numeric(p)
      while(Taux < Trunc-1){
        index <- rcat(prob=probs[1:newvalueindex])
        if(index==newvalueindex){# sample from G_0
          for(j in 1:p){ 
            model$simulate(tildevarNames[j])
            paramaux[j] <- values(model, tildevarNames[j])[1]  
          }
        }else{# sample one of the existing values
          for(j in 1:p){
            paramaux[j] <- uniqueValues[index, j] 
          }
        }
        condaux <- uniqueValues[1:Taux, 1] == paramaux[1] # check if we sample a new atom or an atom that ir in G already
        if(sum(condaux) >0){ # the atom already exists and we have to update the weights and not increase Trunc
          repindex=1
          while(condaux[repindex]==FALSE){
            repindex=repindex+1
          }
          v1prod <- v1prod*(1-vaux)
          vaux <-rbeta(1, 1, conciter+N)
          mv['G', iiter][repindex, 1] <<- mv$G[[iiter]][repindex, 1] + vaux*v1prod
        }else{ # agument the truncation and keep the same parameters
          Taux <-Taux+1
          for(j in 1:p){
            mv$G[[iiter]][Taux, j+1] <<- paramaux[j]
          }
          if( Taux==1 ){
            mv['G', iiter][Taux, 1] <<-vaux
          }else{
            v1prod <- v1prod*(1-vaux)
            vaux <- rbeta(1, 1, conciter+N)
            mv['G', iiter][Taux, 1] <<-vaux*v1prod
          }
        }
      }
      # complete the vector of probabiities and atoms
      mv['G', iiter][Trunc, 1] <<- 1- sum(mv['G', iiter][1:(Trunc-1), 1])
      for(j in 1:p){ 
        model$simulate(tildevarNames[j])
        mv['G', iiter][Trunc, j+1] <<- values(model, tildevarNames[j])[1]  
      }
    }
  },
  
  methods = list( reset = function () {} )
)

