nimMaxLik<- nimbleFunction(
  setup = function(model, paramNodes, calcNodes) {
    ## The following model processing code ensures that any deterministic nodes
    ## between paramNodes and calcNodes are added to calcNodes.
    ## Step 1: Deterministic dependencies of paramNodes
    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
    ## Step 2: Weed out any paramDeps that are already included in calcNodes and/or do not lead to part of calcNodes
    if(length(paramDeps) > 0) {
      keep_paramDeps <- logical(length(paramDeps))
      for(i in seq_along(paramDeps)) {
        if(any(paramDeps[i] == calcNodes)) keep_paramDeps[i] <- FALSE
        else {
          nextDeps <- model$getDependencies(paramDeps[i])
          keep_paramDeps[i] <- any(nextDeps %in% calcNodes)
        }
      }
      paramDeps <- paramDeps[keep_paramDeps]
    }
    ## Add any additional nodes to calcNodes
    calcNodes <- model$expandNodeNames(c(paramDeps, calcNodes), sort = TRUE)
    ## Automated transformation for parameters
    paramsTransformation <- parameterTransform(model, paramNodes)
    pTransform_length <- paramsTransformation$getTransformedLength()
    pTransform_indices <- if(pTransform_length > 1) 1:pTransform_length else c(1, -1)
    one_time_fixes_done <- FALSE
    
    makeUpdateNodes <- nimble:::makeDerivsInfo(wrtNodes = paramNodes, calcNodes = calcNodes, model = model)
    updateNodes     <- makeUpdateNodes$updateNodes
    constantNodes   <- makeUpdateNodes$constantNodes
  },
  run = function(){},
  methods = list(
    fix_one_vec = function(x = integer(1)) {
      if(length(x) == 2) {
        if(x[2] == -1) {
          ans <- integer(length = 1, value = x[1])
          return(ans)
        }
      }
      return(x)
      returnType(integer(1))
    },
    one_time_fixes = function() {
      ## Run this once after compiling; remove extraneous -1 if necessary
      if(one_time_fixes_done) return()
      if(pTransform_length == 1) {
        pTransform_indices <<- fix_one_vec(pTransform_indices)
      }
      one_time_fixes_done <<- TRUE
    },
    ## Log-likelihood 
    logLik = function(pTransform = double(1)) {
      p <- paramsTransformation$inverseTransform(pTransform)
      values(model, paramNodes) <<- p
      ans <- model$calculate(calcNodes) + paramsTransformation$logDetJacobian(pTransform)
      return(ans)
      returnType(double())
    },
    # Gradient of the log-likelihood w.r.t parameters
    gr_logLik_internal = function(pTransform = double(1)) {
      ans <- derivs(logLik(pTransform), wrt = pTransform_indices,
                    order = 1, model = model,
                    updateNodes   = updateNodes, 
                    constantNodes = constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    gr_logLik = function(pTransform = double(1)) {
      ans <- derivs(gr_logLik_internal(pTransform), wrt = pTransform_indices,
                    order = 0, model = model,
                    updateNodes   = updateNodes,
                    constantNodes = updateNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Hessian matrix
    hess_internal = function(pTransform = double(1)) {
      ans <- derivs(gr_logLik_internal(pTransform), wrt = pTransform_indices,
                    order = 1, model = model,
                    updateNodes   = updateNodes, 
                    constantNodes = constantNodes)
      return(ans$jacobian)
      returnType(double(2))
    },
    ## Double tapping
    hess = function(pTransform = double(1)) {
      ans <- derivs(hess_internal(pTransform), wrt = pTransform_indices,
                    order = 0, model = model,
                    updateNodes   = updateNodes,
                    constantNodes = constantNodes)
      hessmat <- matrix(ans$value, nrow = pTransform_length)
      return(hessmat)
      returnType(double(2))
    },
    ## Maximize the log-likelihood
    maxLik = function(pStart = double(1)) {
      pStartTransform <- paramsTransformation$transform(pStart)
      optimControl <- nimOptimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 1000
      optRes <- optim(pStartTransform, logLik, gr_logLik, #hess_internal,
                      method = "BFGS", control = optimControl)
      ## Need to consider the case where optim does not converge
      if(optRes$convergence != 0) {
        print("Warning: optim does not converge properly.")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    }
  ),
  buildDerivs = list(logLik             = list(),
                     gr_logLik_internal = list(),
                     hess_internal      = list())
)
