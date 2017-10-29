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

