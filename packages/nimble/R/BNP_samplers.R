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

# eventhough in the definition of the model N, N2, and N3 can be non admisible, the user can get to this step without realising that.
# the creation of the model sends an ERRROR, but it is builded and then the model can be compiled 
# N: number of observation; N2: number of xi; N3: number of tilde variables
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
    
    # first check that the sampler can be used: N=N2
    N2 <- length(model[[target]]) # number of xi in the model
    N <- length(model$getDependencies(targetElements, dataOnly = TRUE) )
    VarNames <- model$getVarNames()
    data <- model$getDependencies(targetElements[1], dataOnly = TRUE) 
    #i <- 1
    #condaux <- TRUE
    #while(condaux && i<=length(VarNames)){
    #  expandVarNames <- model$expandNodeNames(VarNames[i])
    #  if(expandVarNames[1] == data){
    #    dataVar=VarNames[i]
    #    condaux <- FALSE
    #  }else{
    #    i=i+1
    #  }
    #}
    #N <- length(model[[dataVar]])
    
    if(N != N2){ stop('length of random indexes and observations has to be the same') }
    
    # finding tilde variables:
    tildevarNames=c()
    itildeVar <- 1
    
    Dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in 1:length(Dep)){ 
      Depi <- Dep[i]
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(Depi))
      expr <- parse(text = expr)[[1]]
      if(is.call(expr) && expr[[1]] == '[' && expr[[3]] == targetElements[1]){
        tildevarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
        itildeVar <- itildeVar+1 
      }
    }
    
    # second check that the sampler can be used: N3 < N2. 
    # N3>N2 might be inefficient, but should still run
    # finding length of tilde variables:
    #targetElements[1] <- n
    #out <- try(model$calculate(model$getDependencies(targetElements[1])))
    #try(model$calculate(model$getDependencies('thetatilde[n]')))
    #try(model$calculate('thetatilde[n]'))
    #try(model[['thetatilde[n]']])
    #out=try(model$getDependencies('thetatilde[n]'))
    #if(out=='try-error' && grep("dynamic index out of bounds", out) == 1){stop('stop')}
    N3 <- c()
    for(i in 1: length(tildevarNames)){
      N3[i] <- length(model[[tildevarNames[i]]])
    }
    
    if(sum(N3==N3[1]) != length(N3)){
      stop('tilde variables have different length')
    }
    for(i in 1:length(N3)){
      if(N3[i] < N2){ stop('length of tilde node has to be at least equal to length of random indexes') }   
    }
    
    #tildeVariable <- c()
    #N3 <- c()
    #N3i <- 1
    #for(i in 1:length(VarNames)){
    #  expandVarNames <- model$expandNodeNames(VarNames[i])
    #  VarExpr <- model$getValueExpr(expandVarNames[1])
    #  if(length(VarExpr)==3){ # thetatilde[xi[i]] or dnorm(mean = mu0, sd = lifted_sqrt_oPtau20_cP); si es poisson tb es de largo 3?
    #    if(VarExpr[[1]]=='[' && VarExpr[[3]]=='xi[1]'){ # length of deterministic nodes theta[i] <- thetatilde[xi[i]]
    #      N3[N3i] <- length(model[[VarExpr[[2]]]])
    #      tildeVariable[N3i] <- VarNames[i]
    #      N3i <- N3i+1
    #    }else{  # length of no deterministic node thetatilde[xi[i]]. Ineficient when the mode has deterministic nodes, but works
    #      depVarNames <- model$getDependencies(VarNames[i], self=FALSE)
    #      if( length(depVarNames)>1 ){# avoids empty depVarNames, e.g. when VarNames[i]='y'
    #        if(sum(depVarNames==data)==1){
    #          N3[N3i] <- length(model[[VarNames[i]]])
    #          tildeVariable[N3i] <- VarNames[i]
    #          N3i <- N3i+1
    #        }
    #      } 
    #    }
    #  }
    #}
     
    # other way to do the same
    #if(nInterm==0){ # no intermediate nodes: for a random node with no reparam
    #  stochDep <- model$getDependencies(targetElements[1], self = FALSE) # y[1]
    #  paramExprs <- nimble:::cc_getNodesInExpr(model$getValueExpr(stochDep)) #e.g., 'mean = theta[xi[i]], var = sigma' 
    #  expr <- parse(text = paramExprs)[[1]]
    #  if(is.call(expr) && expr[[1]] == '[') {
    #    N3 <- length(model$expandNodeNames(expr[[2]])) 
    #  } # la  varianza siempre tiene uno deterministico
    #}
    #if(nInterm==1){ # puede ser theta[i], reparam varianza en dnorm(1,s2tilde[xi[i]]). Que pasa si exp(theta[i])?
    #  detDep <- model$getDependencies(targetElements[1], determOnly = TRUE)  # e.g., 'theta[1]'
    #  detDepExpr <- model$getValueExpr(detDep)
    #  if(is.call(detDepExpr) ){ # tested!
    #    if(detDepExpr[[1]] == '['){
    #      stochNodes <- model$expandNodeNames(detDepExpr[[2]])
    #      N3 <- numeric(length(stochNodes))
    #      for(i in 1:length(stochNodes)){
    #        if(is.numeric(model[[stochNodes[i]]])){N3[i] <- 1}
    #      }
    #      N3 <- sum(N3)
    #    }
    #    if(detDepExpr[[1]] == 'sqrt'){ # still have to test it! # case when we have y[i] ~ dnorm(0, var=s2tilde[xi[i]]):  nInterm=1 and detDepExpr = sqrt(s2tilde[xi[1]])
    #      detDepExpr2 <-  detDepExpr[[2]] #detDepExpr2[[1]]=='['????
    #      if(detDepExpr2[[1]] == '['){
    #        stochNodes <- model$expandNodeNames(detDepExpr2[[2]])
    #        N3 <- length(stochNodes) # there can be other 
    #      }
    #    }
    #  }
    #}
    #if(nInterm>=2){
    #  detDep <- model$getDependencies(targetElements[1], determOnly = TRUE)
    #  xiFound <- TRUE
    #  i <- 1
    #  while(xiFound & i<=length(detDep)){
    #    detDepi <- detDep[i]
    #    detDepExpr <- model$getValueExpr(detDepi)
    #    if(is.call(detDepExpr) && detDepExpr[[3]] == 'xi[1]') {
    #      N3 <- length(model$expandNodeNames(detDepExpr[[2]])) 
    #      xiFound <- FALSE
    #    }else{
    #      i <- i+1
    #    }
    #  }
    #}
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
                                        #if(nInterm==0){ # find tilde variables?
        #  intermNodes[i] <- tildeVariable[xi[i]]?????
        #}  
              }
          }
      }
          
    ## determination of conjugacy for one tilde node that has 1 or 0 determnistic nodes
    ## does not find conjugacy when we have random mean and variance defined by deterministic nodes
    ## does find conjugacy when we have random mean  defined or not by deterministic nodes
    ## dont know what does when y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]]). here nInterm=1
    ## dont know what does when y[i] ~ dnorm(0, var=s2tilde[xi[i]]). here nInterm=1
    
    # new checking for conjugacy using tilde variables:
    if(length(tildevarNames)>1){ conjugate <- FALSE }
    if(length(tildevarNames)==1){
      clusterNodes <- model$expandNodeNames(tildevarNames[1])  # e.g., 'thetatilde[1]',...,
      conjugacy <- model$checkConjugacy(clusterNodes[1], restrictLink = 'identity')
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
    cat("Conjugacy detected: ", conjugate)
    
    #if(nInterm <= 1 && type == 'indivCalcs') {  # only detect conjugacy in simple settings for now
    #  conjugate <- TRUE  # at various steps below, may determine this not to be the case
    #  if(nInterm) {  # 1 intermediate; e.g. theta[i] <- thetatilde[xi[i]]
    #    detDep <- model$getDependencies(targetElements[1], determOnly = TRUE)  # e.g., 'theta[1]'
    #    detDepExpr <- model$getValueExpr(detDep)  # e.g., thetatilde[xi[1]]
    #    if(is.call(detDepExpr)){ #detects clusterNodes for thetatilde and s2tilde when nInterm=1
    #      if(detDepExpr[[1]] == '['){  # detects thetatilde when thetatilde defined a deterministc node
    #        if(detDepExpr[[3]] != targetElements[1]) {
    #          conjugate <- FALSE
    #        } else clusterNodes <- model$expandNodeNames(deparse(detDepExpr[[2]]))  # e.g., 'thetatilde[1]',...,
    #      }else{ # detects s2tilde when no deteministic nodes involve s2tilde
    #        expr <- nimble:::cc_getNodesInExpr(detDepExpr) # s2tilde[xi[1]]
    #        if(expr==detDepExpr[[2]]){ # making sure that var= s2tilde[xi[1]] and not var=2+ s2tilde[xi[1]] 
    #          expr <- parse(text = expr)[[1]]
    #          if(is.call(expr) && expr[[1]] == '[' && expr[[3]] != targetElements[1]) {
    #            conjugate <- FALSE
    #          }else{
    #            clusterNodes <- model$expandNodeNames(deparse(expr[[2]]))
    #          }
    #        }else{
    #          conjugate <- FALSE
    #        }
    #      }
    #      #if(detDepExpr[[1]] == 'sqrt'){ # case when we have y[i] ~ dnorm(0, var=s2tilde[xi[i]]):  nInterm=1 and detDepExpr = sqrt(s2tilde[xi[1]]))
    #      #  detDepExpr2 <-  detDepExpr[[2]]
    #      #  if(detDepExpr2[[1]] == '['){
    #      #    if(detDepExpr2[[3]] != targetElements[1]) {
    #      #      conjugate <- FALSE
    #      #    } else clusterNodes <- model$expandNodeNames(deparse(detDepExpr2[[2]]))  # e.g., 'thetatilde[1]',...,
    #      #  }
    #      #}
    #    }else{
    #      conjugate <- FALSE # we should determine conjugacy!?
    #    }
    #  } else {  # 0 intermediates; e.g., y[i] ~ dnorm(theta[xi[i]], 1)
    #    stochDep <- model$getDependencies(targetElements[1], self = FALSE)  # e.g., y[1]
    #    paramExprs <- nimble:::cc_getNodesInExpr(model$getValueExpr(stochDep)) #e.g., 'mean = theta[xi[i]], var = sigma' 
    #    candidates <- lapply(paramExprs, function(expr) {
    #      expr <- parse(text = expr)[[1]]
    #      if(is.call(expr) && expr[[1]] == '[' && expr[[3]] == targetElements[1]) {
    #        model$expandNodeNames(deparse(expr[[2]]))  # e.g., 'theta[1]', ...
    #      } else NULL
    #    })
    #    nonNull <- which(!sapply(candidates, is.null))
    #    if(length(nonNull) != 1) { # 'xi[1]' appears in multiple parameters of 'y[1]'
    #      conjugate <- FALSE
    #    } else clusterNodes <- candidates[[nonNull]]
    #  }
    #  if(conjugate) {  # not ruled out yet
    #    if(length(unique(model$getDeclID(clusterNodes))) != 1) { # cluster nodes are not all from same prior
    #      conjugate <-  FALSE
    #    } else {
    #      conjugacy <- model$checkConjugacy(clusterNodes[1], restrictLink = 'identity')
    #      ## 'restrictLink' is because we don't yet want to deal with cases like:
    #      ## y[i] ~ dnorm(5*theta[i], tau2)
    #      ## as the 5 complicates the marginal
    #      if(length(conjugacy)) {
    #        conjugacy <- conjugacy[[1]]  # it's a one-element list of lists, so simplify it
    #        conjugate <- TRUE
    #        marginalizedParam <- conjugacy$target
    #        ## print(conjugacy) # will show information you can use for conjugate sampler
    #        ## I think you will simply want to use things like this in run code:
    #        ## this would be the case for the simple normal-normal conjugacy
    #        ## getParam(marginalizedParam, 'mean')
    #        ## getParam(marginalizedParam, 'var') + getParam(dataNodes[i], 'var')
    #        ## to construct the marginal you need
    #      } else conjugate <- FALSE
    #    }
    #  }
    #  cat("Conjugacy detected: ", conjugate)
    #} 
    
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
          curLogProb[j] <<- log(conc) + model$getLogProb(dataNodes[i]) #<<-
        }else{
          model[[target]][i] <<- model[[target]][j]
          if(type == 'indivCalcs') {
            if(nInterm >= 1) model$calculate(intermNodes[i])
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



#-- Standarized output nimbleFunction new version
#-- user gives us: 
#   - model: the model
#   - mvSaved: uncompiled modelValues object with the samples of xi, the tilde variables, and the conc parameter, if sampled 
# works fine for modelos with and without deterministic 1 and 2 nodes

# models that work with this sampler:
#-- 1
#   Code=nimbleCode( {
#    for(i in 1:N3){
#      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
#      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0)     }
#    xi[1:N2] ~ dCRP(conc)
#    for(i in 1:N){
#      theta[i] <- thetatilde[xi[i]]
#      s2[i] <- s2tilde[xi[i]]
#      y[i] ~ dnorm(theta[i], var=s2[i])    }
#    conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40  })
#-- 2
#Code=nimbleCode({
#    for(i in 1:N3){
#      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
#      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0)     }
#    xi[1:N2] ~ dCRP(conc)
#    for(i in 1:N){
#      y[i] ~ dnorm(thetatilde[xi[i]], var=s2tilde[xi[i]])#    }
#    conc<-1;mu0<-0; tau20<-40 ; a0<-1; b0<-0.5;   })
#-- 3: as long as xi is monitored!
#   Code=nimbleCode( {
#    for(i in 1:N3){
#      thetatilde[i] ~ dnorm(mean=mu0, var=tau20) 
#      s2tilde[i] ~ dinvgamma(shape=a0, scale=b0)     }
#    xi[1:N2] ~ dCRP(conc)
#    conc ~ dgamma(1, 1)
#    for(i in 1:N){
#      theta[i] <- thetatilde[xi[i]]
#      s2[i] <- s2tilde[xi[i]]
#      y[i] ~ dnorm(theta[i], var=s2[i])    }
#    conc<-1; a0<-1; b0<-0.5; mu0<-0; tau20<-40  })


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
    tildevarNames=c()
    itildeVar <- 1
    stochDepOnly <- FALSE
    
    Dep <- model$getDependencies(targetElements[1], self=FALSE)
    for(i in 1:length(Dep)){ 
      Depi <- Dep[i]
      expr <- nimble:::cc_getNodesInExpr(model$getValueExpr(Depi))
      expr <- parse(text = expr)[[1]]
      if(is.call(expr) && expr[[1]] == '[' && expr[[3]] == targetElements[1]){
        tildevarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
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
    #      tildevarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
    #      itildeVar <- itildeVar+1 
    #    }
        #if(is.call(detDepExpr) && detDepExpr[[1]]=='['){ # deter. nodes that can depend on xi
        #  if(detDepExpr[[3]] == 'xi[1]'){
        #    tildevarNames[itildeVar] <- VarNames[which(VarNames==detDepExpr[[2]])]
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
    #        tildevarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
    #        itildeVar <- itildeVar+1 
    #      }
    #    }
    #  }
    #  for(i in 1:length(paramExprs2)){ # we get the tilde mean
    #    expr <- parse(text = paramExprs2[i])[[1]]
    #    if(is.call(expr) && expr[[1]]=='['){ # case where we have random sigm yand its reparametrication; can it be >3?
    #      if(expr[[3]] == 'xi[1]'){
    #        tildevarNames[itildeVar] <- VarNames[which(VarNames==expr[[2]])]
    #        itildeVar <- itildeVar+1 
    #      }
    #    }
    #  }
    #}
    #tildevarNames <- unique(tildevarNames)
    
    p <- length(tildevarNames)
    
    # storaging object:
    mvConf <- modelValuesConf(vars = 'G', type = 'double', size = list(G = c(Trunc, p+1)) )
    mv <- modelValues(mvConf, m = 1)
    getTildeVarList <- nimbleFunctionList(getTildeVarVirtual)
    #getTildeVarList <- nimble:::nimbleFunctionList(getTildeVarVirtual)
    for(j in 1:p){
      getTildeVarList[[j]] <- getTildeVar(mvSaved, tildevarNames[j])
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
        model$simulate(tildevarNames)
        for(j in 1:p){ 
          paramaux[j] <- values(model, tildevarNames)[j]#values(model, tildevarNames)[(j-1)*N+1]  #
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
          model$simulate(tildevarNames)
          for(j in 1:p){ 
            paramaux[j] <- values(model, tildevarNames)[j]#values(model, tildevarNames)[(j-1)*N+1]  #
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
      model$simulate(tildevarNames)
      for(j in 1:p){ 
        mv['G', iiter][Trunc, j+1] <<- values(model, tildevarNames)[j]#values(model, tildevarNames)[(j-1)*N+1]  
      }
    }
  },
  methods = list( reset = function () {} )
)




#-- Sampler for concentration parameter, conc, of the dCRP distribution.

sampler_Augmented_BetaGamma <- nimbleFunction(
  name = 'sampler_Augmented_BetaGamma',
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



