## NIMBLE Laplace approximation
## Laplace base class
AGHQuad_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
    calcLogLik1 = function(p = double(1)){
      returnType(double())
    },
    calcLogLik2 = function(p = double(1)){
      returnType(double())
    },
    calcLogLik3 = function(p = double(1)){
      returnType(double())
    },
    gr_logLik1 = function(p = double(1)){
      returnType(double(1))
    },
    gr_logLik2 = function(p = double(1)){
      returnType(double(1))
    },
    gr_logLik3 = function(p = double(1)){
      returnType(double(1))
    },
    negHess = function(p = double(1), reTransform = double(1)){
      returnType(double(2))
    },
    update_max_inner_logLik = function(p = double(1)){
      returnType(double(1))
    },
    update_max_inner_logLik_internal = function(p = double(1)){
      returnType(double(1))
    },
    hess_joint_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)){
      returnType(double(2))
    },
    hess_joint_logLik_wrt_p_wrt_re_internal = function(p = double(1), reTransform = double(1)){
      returnType(double(2))
    }
  )
)

## A single Laplace approximation for only one scalar random effect node
buildOneLaplace1D <- function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
  buildOneAGHQuad1D(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart)
}

buildOneAGHQuad1D <- nimbleFunction(
  contains = AGHQuad_BASE,
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
    ## Check the number of random effects is 1
    nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))
    if(length(nre) != 1) stop("Number of random effects for buildOneAGHQuad1D or buildOneLaplace1D must be 1.")
    ## Check and add necessary upstream deterministic nodes into calcNodes
    ## This ensures that deterministic nodes between paramNodes and calcNodes are used.
    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
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
    innerCalcNodes <- calcNodes
    calcNodes <- model$expandNodeNames(c(paramDeps, calcNodes), sort = TRUE)
    wrtNodes <- c(paramNodes, randomEffectsNodes)
    ## Indices of randomEffectsNodes and paramNodes inside wrtNodes
    npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))
    re_indices <- as.numeric(c(npar+1, -1))
    if(npar > 1) p_indices <- as.numeric(1:npar)
    else p_indices <- as.numeric(c(1, -1))
    ## Indices of randomEffectsNodes inside randomEffectsNodes for use in getting the derivative of
    ## the inner log-likelihood (paramNodes fixed) w.r.t. randomEffectsNodes.
    re_indices_inner <- as.numeric(c(1, -1))
    p_and_re_indices <- as.numeric(1:(npar + 1))
    
    ## Set up start values for the inner optimization of Laplace approximation
    if(identical(optimStart, "last")) {
      startID <- 1
      optStart <- numeric(2)
    }
    else if(identical(optimStart, "last.best")) {
      startID <- 2
      optStart <- numeric(2)
    }
    else {
      startID <- 3
      optStart <- as.numeric(c(optimStart, -1))
    }
    
    ## Update and constant nodes for obtaining derivatives using AD
    inner_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes
    
    ## Automated transformation for random effects to ensure range of (-Inf, Inf) 
    reTrans <- parameterTransform(model, randomEffectsNodes)
    
    ## The following are used for caching values and gradient in the Laplace3 system
    logLik3_saved_value <- numeric(1)
    logLik3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    logLik3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    ## The following are used for caching values and gradient in the Laplace3 system
    max_inner_logLik_saved_par <- as.numeric(c(1, -1))
    max_inner_logLik_saved_value <- numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    ## Record the maximum Laplace loglikelihood value for obtaining inner optimization start values
    max_logLik <- -Inf
    max_logLik_saved_re_value <- as.numeric(c(1, -1))
    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE
  },
  run = function(){},
  methods = list(
    fix_one_vec = function(x = double(1)) {
      if(length(x) == 2) {
        if(x[2] == -1) {
          ans <- numeric(length = 1, value = x[1])
          return(ans)
        }
      }
      return(x)
      returnType(double(1))
    },
    set_reInit = function(re = double(1)) {
      reInitTrans <- reTrans$transform(re)
      max_inner_logLik_saved_par <<- reInitTrans
    },
    get_reInitTrans = function() {
      if(startID == 1) ans <- max_inner_logLik_saved_par
      else if(startID == 2) ans <- max_logLik_saved_re_value
      else ans <- reTrans$transform(optStart)
      return(ans)
      returnType(double(1))
    },
    one_time_fixes = function() {
      ## Run this once after compiling; remove extraneous -1 if necessary
      if(one_time_fixes_done) return()
      re_indices <<- fix_one_vec(re_indices)
      re_indices_inner <<- fix_one_vec(re_indices_inner)
      max_inner_logLik_saved_par <<- fix_one_vec(max_inner_logLik_saved_par)
      max_logLik_saved_re_value <<- fix_one_vec(max_logLik_saved_re_value)
      if(startID == 3) optStart <<- fix_one_vec(optStart)
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        logLik3_saved_gr <<- fix_one_vec(logLik3_saved_gr)
        logLik3_previous_p <<- fix_one_vec(logLik3_previous_p)
        max_inner_logLik_previous_p <<- fix_one_vec(max_inner_logLik_previous_p)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit(reInit)
      one_time_fixes_done <<- TRUE
    },
    ## Joint log-likelihood with values of parameters fixed: used only for inner optimization
    inner_logLik = function(reTransform = double(1)) {
      re <- reTrans$inverseTransform(reTransform)
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(innerCalcNodes) + reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    # Gradient of the joint log-likelihood (p fixed) w.r.t. transformed random effects: used only for inner optimization
    gr_inner_logLik_internal = function(reTransform = double(1)) {
      ans <- derivs(inner_logLik(reTransform), wrt = re_indices_inner, order = 1, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    gr_inner_logLik = function(reTransform = double(1)) {
      ans <- derivs(gr_inner_logLik_internal(reTransform), wrt = re_indices_inner, order = 0, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Solve the inner optimization for Laplace approximation
    max_inner_logLik = function(p = double(1)) {
      values(model, paramNodes) <<- p
      model$calculate(paramDeps)
      reInitTrans <- get_reInitTrans()
      fn_init <- inner_logLik(reInitTrans)
      if((fn_init == Inf) | (fn_init == -Inf) | (is.nan(fn_init)) | (is.na(fn_init))) {
        optRes <- optimResultNimbleList$new()
        optRes$par <- reInitTrans
        optRes$value <- -Inf
        optRes$convergence <- -1
        return(optRes)
      }
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik, method = optimMethod, control = optimControl)
      if(optRes$convergence != 0){
        print("Warning: optim does not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Inner optimization using single-taped gradient
    max_inner_logLik_internal = function(p = double(1)) {
      values(model, paramNodes) <<- p
      model$calculate(paramDeps)
      reInitTrans <- get_reInitTrans()
      fn_init <- inner_logLik(reInitTrans)
      if((fn_init == Inf) | (fn_init == -Inf) | (is.nan(fn_init)) | (is.na(fn_init))) {
        optRes <- optimResultNimbleList$new()
        optRes$par <- reInitTrans
        optRes$value <- -Inf
        optRes$convergence <- -1
        return(optRes)
      }
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik_internal, method = optimMethod, control = optimControl)
      if(optRes$convergence != 0){
        print("Warning: optim does not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## These two update methods for max_inner_logLik use the same member data caches
    update_max_inner_logLik = function(p = double(1)) {
      optRes <- max_inner_logLik(p)
      max_inner_logLik_saved_par <<- optRes$par
      max_inner_logLik_saved_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_saved_par)
      returnType(double(1))
    },
    update_max_inner_logLik_internal = function(p = double(1)) {
      optRes <- max_inner_logLik_internal(p)
      max_inner_logLik_saved_par <<- optRes$par
      max_inner_logLik_saved_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_saved_par)
      returnType(double(1))
    },
    ## Joint log-likelihood in terms of parameters and transformed random effects
    joint_logLik = function(p = double(1), reTransform = double(1)) {
      re <- reTrans$inverseTransform(reTransform)
      values(model, paramNodes) <<- p
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(calcNodes) +  reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    ## 1st order partial derivative w.r.t. parameters
    gr_joint_logLik_wrt_p_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik(p, reTransform), wrt = p_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = p_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_updateNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 1st order partial derivative w.r.t. transformed random effects
    gr_joint_logLik_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 2nd order mixed partial derivative w.r.t. parameters and transformed random effects
    hess_joint_logLik_wrt_p_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    hess_joint_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      derivmat <- matrix(value = ans$value, nrow = npar)
      return(derivmat)
      returnType(double(2))
    },
    ## Negative Hessian: 2nd order unmixed partial derivative w.r.t. transformed random effects
    negHess_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(-ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    negHess = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(negHess_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      neghess <- matrix(ans$value, nrow = nre)
      return(neghess)
      returnType(double(2))
    },
    ## Logdet negative Hessian
    logdetNegHess = function(p = double(1), reTransform = double(1)) {
      negHessian <- negHess(p, reTransform)
      ans <- log(negHessian[1,1])
      return(ans)
      returnType(double())
    },
    ## Gradient of logdet (negative) Hessian w.r.t. parameters
    gr_logdetNegHess_wrt_p_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = p_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_p_internal(p, reTransform), wrt = p_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Gradient of logdet (negative) Hessian w.r.t. transformed random effects
    gr_logdetNegHess_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_re_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Put everything (gradient and Hessian) together for Laplace3
    joint_logLik_with_grad_and_hess = function(p = double(1), reTransform = double(1)) {
      # This returns a vector of  concatenated key quantities (see comment below for details)
      # reTransform is the arg max of the inner logLik
      # We could consider returning only upper triangular elements of chol(-Hessian),
      # and re-constituting as a matrix when needed.
      joint_logLik_res <- derivs(joint_logLik(p, reTransform), wrt = p_and_re_indices, order = c(1, 2),
                                 model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      negHessValue <- -joint_logLik_res$hessian[npar + 1, npar + 1, 1]
      logdetNegHessAns <- log(negHessValue)
      hess_wrt_p_wrt_re <- matrix(init = FALSE, nrow = npar, ncol = nre)
      for(i in 1:npar){
        hess_wrt_p_wrt_re[i, 1] <- joint_logLik_res$hessian[i, npar + 1, 1]
      }
      ans <- c(joint_logLik_res$jacobian[1, 1:npar], logdetNegHessAns, negHessValue, hess_wrt_p_wrt_re)
      ## If cholNegHess is considered, indices to components are:
      ## gr_joint_logLik_wrt_p = (1:npar)                    [size = npar]
      ## logdetNegHess         = npar + 1                    [1]
      ## cholNegHess           = npar + 1 + (1 : nre*nre)    [nre x nre]
      ## hess_wrt_p_wrt_re     = npar + 1 + nre*nre + (1:npar*nre)  [npar x nre]
      return(ans)
      returnType(double(1))
    },
    joint_logLik_with_higher_derivs = function(p = double(1), reTransform = double(1)) {
      # value gives results from joint_logLik_with_grad_and_hess
      # jacobian gives derivs of these outputs wrt (p, re).
      # We only need gradient of logdetNegHess, which is the
      #   (1 + npar + 1, given in that order for sanity) row of jacobian
      # Other rows of the jacobian are wasted, but when this function
      # is meta-taped and optimized (part of CppAD), those calculations should be omitted
      higher_order_deriv_res <- derivs(joint_logLik_with_grad_and_hess(p, reTransform), wrt = p_and_re_indices, order = c(0, 1),
                                       model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      ans <- c(higher_order_deriv_res$value, higher_order_deriv_res$jacobian[npar + 1,])
      return(ans)
      returnType(double(1))
    },
    update_logLik3_with_gr = function(p = double(1), reset = logical(0, default = FALSE)) {
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      ans <- derivs(joint_logLik_with_higher_derivs(p, reTransform), wrt = p_and_re_indices, order = 0,
                    model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      ind <- 1
      # all "logLik" here is joint log likelihood (i.e. for p and re)
      gr_logLik_wrt_p <- numeric(value = ans$value[(ind):(ind + npar - 1)], length = npar)
      ind <- ind + npar
      logdetNegHess_value <- ans$value[ind]
      ind <- ind + 1
      # chol_negHess <- matrix(ans$value[(ind):(ind + nre*nre - 1)], nrow = nre, ncol = nre)
      negHessValue <- ans$value[ind]
      ind <- ind + 1
      hess_cross_terms <- numeric(value = ans$value[(ind):(ind + npar*1 - 1)], length = npar*1)
      ind <- ind + npar*1
      gr_logdetNegHess_wrt_p_v <- numeric(value = ans$value[(ind):(ind + npar - 1)], length = npar)
      ind <- ind + npar
      gr_logdetNegHess_wrt_re_v <- ans$value[ind]
      
      logLik_value <- maxValue - 0.5 * logdetNegHess_value + 0.5 * 1 * log(2*pi)
      logLik3_saved_value <<- logLik_value
      
      gr_logLik_v <- gr_logLik_wrt_p - 0.5*(gr_logdetNegHess_wrt_p_v + hess_cross_terms * (gr_logdetNegHess_wrt_re_v / negHessValue))
      logLik3_saved_gr <<- gr_logLik_v
      return(ans$value)
      returnType(double(1))
    },
    logLik3_update = function(p = double(1)) {
      if(any(p != logLik3_previous_p)) {
        update_logLik3_with_gr(p)
        logLik3_previous_p <<- p
      }
    },
    calcLogLik3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      logLik3_update(p)
      if(logLik3_saved_value > max_logLik) {
        max_logLik <<- logLik3_saved_value
        max_logLik_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(logLik3_saved_value)
      returnType(double())
    },
    gr_logLik3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      logLik3_update(p)
      return(logLik3_saved_gr)
      returnType(double(1))
    },
    ## Laplace approximation 2: double tapping with separate components
    calcLogLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ## Laplace approximation
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * 1 * log(2*pi)
      if(ans > max_logLik) {
        max_logLik <<- ans
        max_logLik_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Laplace approximation 1: single tapping with separate components
    calcLogLik1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ## Laplace approximation
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * 1 * log(2*pi)
      if(ans > max_logLik) {
        max_logLik <<- ans
        max_logLik_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation (version 2) w.r.t. parameters
    gr_logLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      negHessian <- negHess(p, reTransform)[1, 1]
      # invNegHessian <- inverse(negHessian)
      grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p(p, reTransform)
      grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re(p, reTransform)[1]
      hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re(p, reTransform)[,1]
      p1 <- gr_joint_logLik_wrt_p(p, reTransform)
      ans <- p1 - 0.5 * (grlogdetNegHesswrtp + hesslogLikwrtpre * (grlogdetNegHesswrtre / negHessian))
      return(ans)
      returnType(double(1))
    },
    ## Gradient of the Laplace approximation (version 1) w.r.t. parameters
    gr_logLik1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_saved_par
      negHessian <- negHess_internal(p, reTransform)[1, 1]
      ## invNegHessian <- inverse(negHessian)
      grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p_internal(p, reTransform)
      grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re_internal(p, reTransform)[1]
      hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform)[,1]
      ans <- gr_joint_logLik_wrt_p_internal(p, reTransform) - 
        0.5 * (grlogdetNegHesswrtp + hesslogLikwrtpre * (grlogdetNegHesswrtre / negHessian))
      return(ans)
      returnType(double(1))
    }
  ),
  buildDerivs = list(inner_logLik                            = list(),
                     joint_logLik                            = list(),
                     gr_joint_logLik_wrt_re                  = list(),
                     negHess                                 = list(),
                     logdetNegHess                           = list(), 
                     gr_inner_logLik_internal                = list(),
                     gr_joint_logLik_wrt_p_internal          = list(),
                     gr_joint_logLik_wrt_re_internal         = list(),
                     hess_joint_logLik_wrt_p_wrt_re_internal = list(),
                     negHess_internal                        = list(),
                     gr_logdetNegHess_wrt_p_internal         = list(),
                     gr_logdetNegHess_wrt_re_internal        = list(),
                     joint_logLik_with_grad_and_hess         = list(ignore = c("i","j")),
                     joint_logLik_with_higher_derivs         = list())
) ## End of buildOneAGHQuad1D


## A single Laplace approximation for models with more than one scalar random effect node
buildOneLaplace <- function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
  buildOneAGHQuad(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart)
}

buildOneAGHQuad <- nimbleFunction(
  contains = AGHQuad_BASE,
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
    ## Check and add necessary (upstream) deterministic nodes into calcNodes
    ## This ensures that deterministic nodes between paramNodes and calcNodes are used.
    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
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
    innerCalcNodes <- calcNodes
    calcNodes <- model$expandNodeNames(c(paramDeps, calcNodes), sort = TRUE)
    wrtNodes <- c(paramNodes, randomEffectsNodes)
    ## Indices of randomEffectsNodes and paramNodes inside wrtNodes
    reTrans <- parameterTransform(model, randomEffectsNodes)
    npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))
    nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))
    nreTrans <- reTrans$getTransformedLength()
    if(nreTrans > 1) reTrans_indices <- as.numeric((npar+1):(npar+nreTrans))
    else reTrans_indices <- as.numeric(c(npar+1, -1)) 
    if(npar > 1) p_indices <- as.numeric(1:npar)
    else p_indices <- as.numeric(c(1, -1))
    ## Indices of randomEffectsNodes inside randomEffectsNodes for use in getting the derivative of
    ## the inner log-likelihood (paramNodes fixed) w.r.t. randomEffectsNodes.
    if(nreTrans > 1) reTrans_indices_inner <- as.numeric(1:nreTrans)
    else reTrans_indices_inner <- as.numeric(c(1, -1))
    p_and_reTrans_indices <- as.numeric(1:(npar + nreTrans))
    
    ## Set up start values for the inner optimization of Laplace approximation
    if(identical(optimStart, "last")) {
      startID <- 1
      optStart <- numeric(nre)
    }
    else if(identical(optimStart, "last.best")) {
      startID <- 2
      optStart <- numeric(nre)
    }
    else {
      startID <- 3
      optStart <- optimStart
    }
    ## Update and constant nodes info for obtaining derivatives using AD
    inner_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes
    
    ## The following are used for caching values and gradient in the Laplace3 system
    logLik3_saved_value <- -Inf #numeric(1)
    logLik3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    logLik3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    
    max_inner_logLik_saved_par <- if(nreTrans > 1) numeric(nreTrans) else as.numeric(c(1, -1))
    max_inner_logLik_saved_value <- -Inf #numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    
    ## Record the maximum Laplace loglikelihood value for obtaining inner optimization start values
    max_logLik <- -Inf
    max_logLik_saved_re_value <- if(nreTrans > 1) numeric(nreTrans) else as.numeric(c(1, -1))
    
    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE
    update_once <- TRUE
    gr_inner_update_once <- TRUE
    gr_inner_logLik_force_update <- TRUE
    gr_inner_logLik_first <- TRUE
    negHess_inner_update_once <- TRUE
    negHess_inner_logLik_force_update <- TRUE
    negHess_inner_logLik_first <- TRUE
  },
  run = function(){},
  methods = list(
    fix_one_vec = function(x = double(1)) {
      if(length(x) == 2) {
        if(x[2] == -1) {
          ans <- numeric(length = 1, value = x[1])
          return(ans)
        }
      }
      return(x)
      returnType(double(1))
    },
    set_reInit = function(re = double(1)) {
      reInitTrans <- reTrans$transform(re)
      max_inner_logLik_saved_par <<- reInitTrans
    },
    get_reInitTrans = function() {
      if(startID == 1) ans <- max_inner_logLik_saved_par
      else if(startID == 2) ans <- max_logLik_saved_re_value
      else ans <- reTrans$transform(optStart)
      return(ans)
      returnType(double(1))
    },
    set_gr_inner_update = function(update = logical(0, default = TRUE)) {
      gr_inner_update_once <<- update
    },
    set_negHess_inner_update = function(update = logical(0, default = TRUE)) {
      negHess_inner_update_once <<- update
    },
    set_params = function(p = double(1)) {
      values(model, paramNodes) <<- p
      model$calculate(paramDeps)
      gr_inner_update_once <<- TRUE
      negHess_inner_update_once <<- TRUE
    },
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(nre == 1) {
        reTrans_indices <<- fix_one_vec(reTrans_indices)
        reTrans_indices_inner <<- fix_one_vec(reTrans_indices_inner)
        max_inner_logLik_saved_par <<- fix_one_vec(max_inner_logLik_saved_par)
        max_logLik_saved_re_value <<- fix_one_vec(max_logLik_saved_re_value)
      }
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        logLik3_saved_gr <<- fix_one_vec(logLik3_saved_gr)
        logLik3_previous_p <<- fix_one_vec(logLik3_previous_p)
        max_inner_logLik_previous_p <<- fix_one_vec(max_inner_logLik_previous_p)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit(reInit)
      one_time_fixes_done <<- TRUE
    },
    ## Joint log-likelihood with values of parameters fixed: used only for inner optimization
    inner_logLik = function(reTransform = double(1)) {
      re <- reTrans$inverseTransform(reTransform)
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(innerCalcNodes) + reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    # Gradient of the joint log-likelihood (p fixed) w.r.t. transformed random effects: used only for inner optimization
    gr_inner_logLik_internal = function(reTransform = double(1)) {
      ans <- derivs(inner_logLik(reTransform), wrt = reTrans_indices_inner, order = 1, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    gr_inner_logLik = function(reTransform = double(1)) {
      ans <- derivs(gr_inner_logLik_internal(reTransform), wrt = reTrans_indices_inner, order = 0, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = gr_inner_logLik_force_update | gr_inner_update_once)
      gr_inner_update_once <<- FALSE
      return(ans$value)
      returnType(double(1))
    },
    negHess_inner_logLik_internal = function(reTransform = double(1)) {
      ans <- derivs(gr_inner_logLik_internal(reTransform), wrt = reTrans_indices_inner, order = 1, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes)
      return(-ans$jacobian)
      returnType(double(2))
    },
    # We also tried double-taping straight to second order. That was a bit slower.
    negHess_inner_logLik = function(reTransform = double(1)) {
      ans <- derivs(negHess_inner_logLik_internal(reTransform), wrt = reTrans_indices_inner, order = 0, model = model,
                    updateNodes = inner_updateNodes, constantNodes = inner_constantNodes,
                    do_update = negHess_inner_logLik_force_update | negHess_inner_update_once)
      negHess_inner_update_once <<- FALSE
      neghess <- matrix(ans$value, nrow = nreTrans)
      return(neghess)
      returnType(double(2))
    },
    record_negHess_inner_logLik = function(reTransform = double(1)) {
      negHess_inner_logLik_force_update <<- TRUE
      negHess_inner_logLik(reTransform) # record
      negHess_inner_logLik_first <<- FALSE
      negHess_inner_logLik_force_update <<- FALSE
    },
    ## Solve the inner optimization for Laplace approximation
    max_inner_logLik = function(p = double(1)) {
      set_params(p)
      reInitTrans <- get_reInitTrans()
      fn_init <- inner_logLik(reInitTrans)
      if((fn_init == Inf) | (fn_init == -Inf) | (is.nan(fn_init)) | (is.na(fn_init))) {
        optRes <- optimResultNimbleList$new()
        optRes$par <- reInitTrans
        optRes$value <- -Inf
        optRes$convergence <- -1
        return(optRes)
      }
      if(gr_inner_logLik_first) { 
        gr_inner_logLik_force_update <<- TRUE
        gr_inner_logLik(reInitTrans) 
        gr_inner_logLik_first <<- FALSE
        gr_inner_logLik_force_update <<- FALSE
      }
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik, method = optimMethod, control = optimControl)
      if(optRes$convergence != 0){
        print("Warning: optim does not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    },
    max_inner_logLik_internal = function(p = double(1)) {
      set_params(p)
      reInitTrans <- get_reInitTrans()
      fn_init <- inner_logLik(reInitTrans)
      if((fn_init == Inf) | (fn_init == -Inf) | (is.nan(fn_init)) | (is.na(fn_init))) {
        optRes <- optimResultNimbleList$new()
        optRes$par <- reInitTrans
        optRes$value <- -Inf
        optRes$convergence <- -1
        return(optRes)
      }
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik_internal, method = optimMethod, control = optimControl)
      if(optRes$convergence != 0){
        print("Warning: optim does not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## These two update methods for max_inner_logLik use the same member data caches
    update_max_inner_logLik = function(p = double(1)) {
      optRes <- max_inner_logLik(p)
      max_inner_logLik_saved_par <<- optRes$par
      max_inner_logLik_saved_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_saved_par)
      returnType(double(1))
    },
    update_max_inner_logLik_internal = function(p = double(1)) {
      optRes <- max_inner_logLik_internal(p)
      max_inner_logLik_saved_par <<- optRes$par
      max_inner_logLik_saved_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_saved_par)
      returnType(double(1))
    },
    ## Joint log-likelihood in terms of parameters and transformed random effects
    joint_logLik = function(p = double(1), reTransform = double(1)) {
      re <- reTrans$inverseTransform(reTransform)
      values(model, paramNodes) <<- p
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(calcNodes) +  reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    ## 1st order partial derivative w.r.t. parameters
    gr_joint_logLik_wrt_p_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik(p, reTransform), wrt = p_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = p_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_updateNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 1st order partial derivative w.r.t. transformed random effects
    gr_joint_logLik_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik(p, reTransform), wrt = reTrans_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = reTrans_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 2nd order mixed partial derivative w.r.t. parameters and transformed random effects
    hess_joint_logLik_wrt_p_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = reTrans_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    hess_joint_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform), wrt = reTrans_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      derivmat <- matrix(value = ans$value, nrow = npar)
      return(derivmat)
      returnType(double(2))
    },
    ## Negative Hessian: 2nd order unmixed partial derivative w.r.t. transformed random effects
    negHess_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = reTrans_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(-ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    negHess = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(negHess_internal(p, reTransform), wrt = reTrans_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes, do_update = update_once)
      # update_once <<- FALSE
      neghess <- matrix(ans$value, nrow = nreTrans)
      return(neghess)
      returnType(double(2))
    },
    reset_update = function(update = logical(0, default = TRUE)) {
      update_once <<- update
    },
    ## Logdet negative Hessian
    logdetNegHess = function(p = double(1), reTransform = double(1)) {
      negHessian <- negHess(p, reTransform)
      cholNegHess <- chol(negHessian)
      ans <- 2 * sum(log(diag(cholNegHess)))
      return(ans)
      returnType(double())
    },
    ## Gradient of logdet (negative) Hessian w.r.t. parameters
    gr_logdetNegHess_wrt_p_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = p_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_p_internal(p, reTransform), wrt = p_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Gradient of logdet (negative) Hessian w.r.t. transformed random effects
    gr_logdetNegHess_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = reTrans_indices, order = 1, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_re_internal(p, reTransform), wrt = reTrans_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Put everything (gradient and Hessian) together for Laplace3
    joint_logLik_with_grad_and_hess = function(p = double(1), reTransform = double(1)) {
      # This returns a vector of  concatenated key quantities (see comment below for details)
      # reTransform is the arg max of the inner logLik
      # We could consider returning only upper triangular elements of chol(-Hessian),
      #  and re-constituting as a matrix when needed.
      joint_logLik_res <- derivs(joint_logLik(p, reTransform), wrt = p_and_reTrans_indices, order = c(1, 2),
                                 model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      negHessUpper <- matrix(init = FALSE, nrow = nre, ncol = nreTrans)
      for(i in 1:nreTrans){
        for(j in i:nreTrans){
          negHessUpper[i,j] <- -joint_logLik_res$hessian[npar + i, npar + j, 1]
        }
      }
      cholNegHess <- chol(negHessUpper)
      logdetNegHessAns <- 2 * sum(log(diag(cholNegHess)))
      hess_wrt_p_wrt_re <- matrix(init = FALSE, nrow = npar, ncol = nre)
      for(i in 1:npar){
        for(j in 1:nreTrans){
          hess_wrt_p_wrt_re[i, j] <- joint_logLik_res$hessian[i, npar + j, 1]
        }
      }
      ans <- c(joint_logLik_res$jacobian[1, 1:npar], logdetNegHessAns, cholNegHess, hess_wrt_p_wrt_re)
      ## Indices to components of this are:
      ## gr_joint_logLik_wrt_p = (1:npar)                    [size = npar]
      ## logdetNegHess         = npar + 1                    [1]
      ## cholNegHess           = npar + 1 + (1 : nreTrans * nreTrans)    [nreTrans x nreTrans]
      ## hess_wrt_p_wrt_re     = npar + 1 + nre*nre + (1:npar*nreTrans)  [npar x nreTrans]
      return(ans)
      returnType(double(1))
      # return a concatenated vector
    },
    joint_logLik_with_higher_derivs = function(p = double(1), reTransform = double(1)) {
      higher_order_deriv_res <- derivs(joint_logLik_with_grad_and_hess(p, reTransform), wrt = p_and_reTrans_indices, 
                                       order = c(0, 1), model = model,
                                       updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      # value gives results from joint_logLik_with_grad_and_hess
      # jacobian gives derivs of these outputs wrt (p, re).
      # We only need gradient of logdetNegHess, which is the
      #   (1 + npar + 1, given in that order for sanity) row of jacobian
      # Other rows of the jacobian are wasted, but when this function
      # is meta-taped and optimized (part of CppAD), those calculations should be omitted
      ans <- c(higher_order_deriv_res$value, higher_order_deriv_res$jacobian[npar + 1,])
      return(ans)
      returnType(double(1))
    },
    update_logLik3_with_gr = function(p = double(1), reset = logical(0, default = FALSE)) {
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      ans <- derivs(joint_logLik_with_higher_derivs(p, reTransform), wrt = p_and_reTrans_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      ind <- 1
      # all "logLik" here is joint log likelihood (i.e. for p and re)
      gr_logLik_wrt_p <- ans$value[(ind):(ind + npar - 1)]
      ind <- ind + npar
      logdetNegHess_value <- ans$value[ind]
      ind <- ind + 1
      chol_negHess <- matrix(ans$value[(ind):(ind + nreTrans*nreTrans - 1)], nrow = nreTrans, ncol = nreTrans)
      ind <- ind + nreTrans*nreTrans
      hess_cross_terms <- matrix(ans$value[(ind):(ind + npar*nreTrans - 1)], nrow = npar, ncol = nreTrans)
      ind <- ind + npar*nreTrans
      gr_logdetNegHess_wrt_p_v <- ans$value[(ind):(ind + npar - 1)]
      ind <- ind + npar
      gr_logdetNegHess_wrt_re_v <- ans$value[(ind):(ind + nreTrans - 1)]
      
      logLik_value <- maxValue - 0.5 * logdetNegHess_value + 0.5 * nreTrans * log(2*pi)
      logLik3_saved_value <<- logLik_value
      
      # We need A^T inverse(negHess) B
      # where A = gr_logdetNegHess_wrt_re_v (a vector treated by default as a one-column matrix)
      #  and  B = t(hess_cross_terms)
      # We avoid forming the matrix inverse because we have negHess = U^T U, where U = chol(negHess)
      #    so inverse(negNess) = inverse(U) inverse(U^T), and inverse(U^T) = inverse(U)^T
      # Since U it upper triangular, it is typically more efficient to do forwardsolve and/or backsolve
      #    than to actually form inverse(U) or inverse(negHess)
      # We have (A^T inverse(U) ) ( inverse(U^T) B) = v^T w
      #    v^T = A^T inverse(U), so v = inverse(U^T) A = fowardsolve(U^T, gr_logdetNegHess_wrt_re_v )
      #    w = inverse(U^T) B, so w = forwardsolve(U^T, t(hess_cross_terms))
      #
      # We could return the chol and hess_cross_terms from the derivs steps
      # in transposed form since that's how we need them here.
      v <- forwardsolve(t(chol_negHess), gr_logdetNegHess_wrt_re_v)
      w <- forwardsolve(t(chol_negHess), t(hess_cross_terms))
      gr_logLik_v <- gr_logLik_wrt_p - 0.5*(gr_logdetNegHess_wrt_p_v + v %*% w )
      # print( gr_logLik_v )
      logLik3_saved_gr <<- numeric(gr_logLik_v, length = npar)
      return(ans$value)
      returnType(double(1))
    },
    logLik3_update = function(p = double(1)) {
      if(any(p != logLik3_previous_p)) {
        update_logLik3_with_gr(p)
        logLik3_previous_p <<- p
      }
    },
    calcLogLik3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      logLik3_update(p)
      ans <- logLik3_saved_value
      if(ans > max_logLik) {
        max_logLik <<- ans
        max_logLik_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    gr_logLik3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      logLik3_update(p)
      return(logLik3_saved_gr)
      returnType(double(1))
    },
    ## Laplace approximation 2: double tapping with separate components
    calcLogLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      if(ans > max_logLik) {
        max_logLik <<- ans
        max_logLik_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Laplace approximation 1: single tapping with separate components
    calcLogLik1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      if(ans > max_logLik) {
        max_logLik <<- ans
        max_logLik_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation 2 w.r.t. parameters
    gr_logLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      negHessian <- negHess(p, reTransform)
      invNegHessian <- inverse(negHessian)
      grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p(p, reTransform)
      grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re(p, reTransform)
      hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re(p, reTransform)
      ans <- gr_joint_logLik_wrt_p(p, reTransform) - 
        0.5 * (grlogdetNegHesswrtp + (grlogdetNegHesswrtre %*% invNegHessian) %*% t(hesslogLikwrtpre))
      return(ans[1,])
      returnType(double(1))
    },
    ## Gradient of the Laplace approximation 1 w.r.t. parameters
    gr_logLik1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_saved_par
      negHessian <- negHess_internal(p, reTransform)
      invNegHessian <- inverse(negHessian)
      grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p_internal(p, reTransform)
      grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re_internal(p, reTransform)
      hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform)
      ans <- gr_joint_logLik_wrt_p_internal(p, reTransform) -
        0.5 * (grlogdetNegHesswrtp + (grlogdetNegHesswrtre %*% invNegHessian) %*% t(hesslogLikwrtpre))
      return(ans[1,])
      returnType(double(1))
    }
  ),
  buildDerivs = list(inner_logLik                            = list(),
                     joint_logLik                            = list(),
                     gr_joint_logLik_wrt_re                  = list(),
                     negHess                                 = list(),
                     logdetNegHess                           = list(), 
                     gr_inner_logLik_internal                = list(),
                     gr_joint_logLik_wrt_p_internal          = list(),
                     gr_joint_logLik_wrt_re_internal         = list(),
                     hess_joint_logLik_wrt_p_wrt_re_internal = list(),
                     negHess_internal                        = list(),
                     gr_logdetNegHess_wrt_p_internal         = list(),
                     gr_logdetNegHess_wrt_re_internal        = list(),
                     joint_logLik_with_grad_and_hess         = list(ignore = c("i","j")),
                     joint_logLik_with_higher_derivs         = list(),
                     negHess_inner_logLik_internal           = list())
) ## End of buildOneAGHQuad

#' Process model to organize nodes for marginalization as by Laplace or adaptive
#' Gauss-Hermite quadrature (AGHQuad) approximations
#' @export
setupMargNodes <- function(model, paramNodes, randomEffectsNodes, calcNodes,
                           calcNodesOther,
                           split = TRUE,
                           warn = TRUE) {
  paramProvided         <- !missing(paramNodes)
  reProvided            <- !missing(randomEffectsNodes)
  calcProvided          <- !missing(calcNodes)
  calcOtherProvided <- !missing(calcNodesOther)
  # We considered a feature to allow params to be nodes without priors. This is a placeholder in case
  # we ever pursue that again.
  allowNonPriors <- FALSE
  # We may need to use determ and stochastic dependencies of parameters multiple times below
  # Define these to avoid repeated computation
  # A note for future: determ nodes between parameters and calcNodes are needed inside buildOneAGHQuad
  # and buildOneAGHQuad1D. In the future, these could be all done here to be more efficient
  paramDetermDeps <- character(0)
  paramStochDeps  <- character(0)
  paramDetermDepsCalculated <- FALSE
  paramStochDepsCalculated  <- FALSE
  
  # 1. Default parameters are stochastic top-level nodes.
  #    We considered an argument allowNonPriors, defaulting to FALSE. If TRUE, the default params would be all top-level stochastic
  #    nodes with no RHSonly nodes as parents and RHSonly nodes (handling of constants TBD, since non-scalars would be converted to data) that have stochastic dependencies
  #    (And then top-level stochastic nodes with RHSonly nodes as parents are essentially latent/data nodes,
  #    some of which would need to be added to randomEffectsNodes below.)
  stochTopNodes2reNodes <- character(0)
  if(!paramProvided) {
    paramNodes <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
  } else {
    paramNodes <- model$expandNodeNames(paramNodes)
  }

  # 2. Default random effects are latent nodes that are stochastic dependencies of params.
  if((!reProvided) || warn) {
    paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
    paramStochDepsCalculated <- TRUE
    latentNodes <- model$getNodeNames(latentOnly = TRUE)
    reNodesDefault <- intersect(latentNodes, paramStochDeps)
  }
  if(reProvided){
    if(is.null(randomEffectsNodes) || isFALSE(randomEffectsNodes)) randomEffectsNodes <- character(0)
    randomEffectsNodes <- model$expandNodeNames(randomEffectsNodes)
  }
  # 3. Optionally check random effects if they were provided (not default)
  if(reProvided && warn) {
    # First check is for random effects that should have been included but weren't
    reCheck <- setdiff(reNodesDefault, randomEffectsNodes)
    if(length(reCheck)) {
      errorNodes <- paste0(head(reCheck, n = 4), sep = "", collapse = ", ")
      if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
      warning(paste0("There are some random effects (latent states) in the model that look\n",
                     "like they should be included in randomEffectsNodes for Laplace or AGHQuad approximation\n",
                     "for the provided (or default) paramNodes:\n",
                     errorNodes, "\n",
                     "To silence this warning, include \'warn = FALSE\' in the control list."))
    }
    # Second check is for random effects that were included but look unnecessary
    reCheck <- setdiff(randomEffectsNodes, reNodesDefault)
    if(length(reCheck)) {
      errorNodes <- paste0(head(reCheck, n = 4), sep = "", collapse = ", ")
      if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
      warning(paste0("There are some randomEffectsNodes provided that look like\n",
                     "they are not needed for Laplace or AGHQuad approximation for the\n",
                     "provided (or default) paramNodes:\n",
                     errorNodes, "\n",
                     "To silence this warning, include \'warn = FALSE\' in the control list."))
    }
  }
  # End of step 2
  if(!reProvided) {
    randomEffectsNodes <- reNodesDefault
  }
  # 4. Default calc nodes are dependencies of random effects
  #    Note that the internal Laplace nimbleFunctions (one for each conditionally independent set)
  #    will fill in deterministic nodes between paramNodes and randomEffectsNodes
  if((!calcProvided) || warn) {
    calcNodesDefault <- model$getDependencies(randomEffectsNodes)
  }
  if(calcProvided){
    if(is.null(calcNodes) || isFALSE(calcNodes)) calcNodes <- character(0)
    calcNodes <- model$expandNodeNames(calcNodes)
  }
  # 5. Optionally check calcNodes if they were provided (not default)
  if(calcProvided && warn) {
    # First check is for calcNodes that look necessary but were omitted
    calcCheck <- setdiff(calcNodesDefault, calcNodes)
    if(length(calcCheck)) {
      errorNodes <- paste0(head(calcCheck, n = 4), sep = "", collapse = ", ")
      if(length(calcCheck) > 4) errorNodes <- paste(errorNodes, "...")
      warning(paste0("There are some model nodes that look like they should be\n",
                     "included in the calcNodes for Laplace or AGHQuad approximation because\n",
                     "they are dependencies of some randomEffectsNodes:\n",
                     errorNodes, "\n",
                     "To silence this warning, include \'warn = FALSE\' in the control list."))
    }
    # Second check is for calcNodes that look unnecessary
    # If some determ nodes between paramNodes and randomEffectsNodes are provided in calcNodes 
    # then that's ok and we should not throw a warning message. 
    calcCheck <- setdiff(calcNodes, calcNodesDefault)
    errorNodes <- calcCheck[model$getNodeType(calcCheck)=="stoch"]
    determCalcCheck <- setdiff(calcCheck, errorNodes)
    numDetermCalcCheck <- length(determCalcCheck)
    # Check other determ nodes
    if(numDetermCalcCheck){
      paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      paramDetermDepsCalculated <- TRUE
      for(i in 1:numDetermCalcCheck){
        if(!(determCalcCheck[i] %in% paramDetermDeps) || 
           !(any(model$getDependencies(determCalcCheck[i], self = FALSE) %in% calcNodesDefault))){
          errorNodes <- c(errorNodes, determCalcCheck[i])
        }
      }
    }
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      warning(paste0("There are some calcNodes provided that look like\n",
                     "they are not needed for Laplace or AGHQuad approximation over\n",
                     "the provided (or default) randomEffectsNodes:\n",
                     outErrorNodes, "\n",
                     "To silence this warning, include \'warn = FALSE\' in the control list."))
    }
  }
  # Finish step 4
  if(!calcProvided){
    calcNodes <- calcNodesDefault
  }
  # 6. Default calcNodesNoLaplce: nodes needed for full model likelihood but
  #    that are not involved in the marginalization done by Laplace.
  #    Default is a bit complicated: All dependencies from paramNodes to
  #    stochastic nodes that are not part of calcNodes. Note that calcNodes
  #    does not necessarily contain deterministic nodes between paramNodes and
  #    randomEffectsNodes. We don't want to include those in calcNodesOther.
  #    (A deterministic that is needed for both calcNodes and calcNodesOther should be included.)
  #    So we have to first do a setdiff on stochastic nodes and then fill in the
  #    deterministics that are needed.
  if(!calcOtherProvided || warn) {
    if(!paramStochDepsCalculated){
      paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
      paramStochDepsCalculated <- TRUE
    }
    calcNodesOtherDefault <- setdiff(paramStochDeps, calcNodes)
  }
  if(calcOtherProvided){
    if(is.null(calcNodesOther) || isFALSE(calcNodesOther)) calcNodesOther <- character(0)
    calcNodesOther <- model$expandNodeNames(calcNodesOther, sort = TRUE)
    if((length(calcNodesOther) > 0) && !any(model$getNodeType(calcNodesOther)=="stoch")){
      warning("There are no stochastic nodes in the calcNodesOther provided for Laplace or AGHQuad approximation.")
    }
  }
  if(!calcOtherProvided){
    calcNodesOther <- calcNodesOtherDefault
  }
  if(calcOtherProvided && warn) {
    calcNoLaplaceCheck <- setdiff(calcNodesOtherDefault, calcNodesOther)
    if(length(calcNoLaplaceCheck)) {
      # We only check missing stochastic nodes; determ nodes will be added below
      missingStochNodesInds <- which((model$getNodeType(calcNoLaplaceCheck)) == "stoch")
      numMissingStochNodes <- length(missingStochNodesInds)
      missingStochNodes <- calcNoLaplaceCheck[missingStochNodesInds]
      if(numMissingStochNodes){
        errorNodes <- paste0(head(missingStochNodes, n = 4), sep = "", collapse = ", ")
        if(numMissingStochNodes > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some model nodes (stochastic) that look like they should be\n",
                       "included in the calcNodesOther for parts of the likelihood calculation\n",
                       "outside of Laplace or AGHQuad approximation:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'warn = FALSE\' in the control list."))
      }
    }
    # Check redundant stochastic nodes
    calcNoLaplaceCheck <- setdiff(calcNodesOther, calcNodesOtherDefault)
    stochCalcNoLaplaceCheck <- calcNoLaplaceCheck[model$getNodeType(calcNoLaplaceCheck)=="stoch"]
    # Check redundant determ nodes
    determCalcNoLaplaceCheck <- setdiff(calcNoLaplaceCheck, stochCalcNoLaplaceCheck)
    numDetermCalcNoLaplaceCheck <- length(determCalcNoLaplaceCheck)
    errorNodes <- character(0)
    if(numDetermCalcNoLaplaceCheck){
      if(!paramDetermDepsCalculated) {
        paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
        paramDetermDepsCalculated <- TRUE
      }
      for(i in 1:numDetermCalcNoLaplaceCheck){
        if(!(determCalcNoLaplaceCheck[i] %in% paramDetermDeps) || 
           !(any(model$getDependencies(determCalcNoLaplaceCheck[i], self = FALSE) %in% calcNodesOtherDefault))){
          errorNodes <- c(errorNodes, determCalcNoLaplaceCheck[i])
        }
      }
    }
    errorNodes <- c(stochCalcNoLaplaceCheck, errorNodes)
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      warning(paste0("There are some nodes provided in calcNodesOther that look like\n",
                     "they are not needed for parts of the likelihood calculation\n",
                     "outside of Laplace or AGHQuad approximation:\n",
                     outErrorNodes, "\n",
                     "To silence this warning, include \'warn = FALSE\' in the control list."))
    }
  }
  # Check and add necessary (upstream) deterministic nodes into calcNodesOther
  # This ensures that deterministic nodes between paramNodes and calcNodesOther are used.
  num_calcNodesOther <- length(calcNodesOther)
  if(num_calcNodesOther > 0){
    if(!paramDetermDepsCalculated) {
      paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      paramDetermDepsCalculated <- TRUE
    }
    numParamDetermDeps <- length(paramDetermDeps)
    if(numParamDetermDeps > 0) {
      keep_paramDetermDeps <- logical(numParamDetermDeps)
      for(i in seq_along(paramDetermDeps)) {
        nextDeps <- model$getDependencies(paramDetermDeps[i])
        keep_paramDetermDeps[i] <- any(nextDeps %in% calcNodesOther)
      }
      paramDetermDeps <- paramDetermDeps[keep_paramDetermDeps]
    }
    calcNodesOther <- model$expandNodeNames(c(paramDetermDeps, calcNodesOther), sort = TRUE)
  }

  # 7. Do the splitting into sets (if given) or conditionally independent sets (if TRUE)
  givenNodes <- NULL
  reSets <- list()
  if(length(randomEffectsNodes)) {
    if(isFALSE(split)) {
      reSets <- list(randomEffectsNodes)
    } else {
      if(isTRUE(split)) {
        givenNodes <- setdiff(c(paramNodes, calcNodes),
                              c(randomEffectsNodes,
                                model$getDependencies(randomEffectsNodes, determOnly=TRUE)))
        reSets <- model$getConditionallyIndependentSets(
          nodes = randomEffectsNodes, givenNodes = givenNodes,
          unknownAsGiven = allowNonPriors)
      }
      else if(is.numeric(split)){
        reSets <- split(randomEffectsNodes, split)
      }
      else stop("Invalid value for \'split\'.")
    }
  }
  list(paramNodes = paramNodes,
       randomEffectsNodes = randomEffectsNodes,
       calcNodes = calcNodes,
       calcNodesOther = calcNodesOther,
       givenNodes = givenNodes,
       randomEffectsSets = reSets
       )
}

## Main function for Laplace approximation
#' @rdname laplace 
#' @export
buildLaplace <- function(model, paramNodes, randomEffectsNodes, calcNodes, calcNodesOther,
                               control = list()) {
 buildAGHQuad(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes, calcNodesOther,
   control)
}

## Main function for Adaptive Gauss-Hermite Quadratrue
buildAGHQuad <- nimbleFunction(
  name = 'AGHQuad',
  setup = function(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes, calcNodesOther,
                   control = list()) {
    if(is.null(control$split)) split <- TRUE else split <- control$split
    if(is.null(control$warn))   warn <- TRUE else  warn <- control$warn
    # Possible future feature
    #if(is.null(control$allowNonPriors)) allowNonPriors <- FALSE else  allowNonPriors <- control$allowNonPriors
    allowNonPriors <- FALSE

    MargNodes <- NULL
    if(!missing(paramNodes)) {
      if(is.list(paramNodes)) {
        # The user called setupMargNodes and provided a list of that format to paramNodes.
        MargNodes <- paramNodes
      }
    }
    if(is.null(MargNodes)) {
      MargNodes <- setupMargNodes(model = model, paramNodes = paramNodes,
                                  randomEffectsNodes = randomEffectsNodes,
                                  calcNodes = calcNodes,
                                  calcNodesOther = calcNodesOther,
                                  # allowNonPriors = allowNonPriors,
                                  split = split,
                                  warn = warn)
    }
    paramNodes <- MargNodes$paramNodes
    randomEffectsNodes <- MargNodes$randomEffectsNodes
    calcNodes <- MargNodes$calcNodes
    calcNodesOther <- MargNodes$calcNodesOther
    num_calcNodesOther <- length(calcNodesOther)
    # MargNodes$randomEffectsSets will be extracted below if needed

    if(length(calcNodesOther)) {
      otherLogLik_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = paramNodes, calcNodes = calcNodesOther)
      otherLogLik_updateNodes   <- otherLogLik_derivsInfo$updateNodes
      otherLogLik_constantNodes <- otherLogLik_derivsInfo$constantNodes
    }
    else { ## calcNodesOther is empty
      otherLogLik_updateNodes   <- character(0)
      otherLogLik_constantNodes <- character(0)
    }
    ## Out and inner optimization settings
    outOptControl   <- nimOptimDefaultControl()
    innerOptControl <- nimOptimDefaultControl()
    optimControlArgNames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", 
                              "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    if(!is.null(control$outOptimControl)){
      validNames <- intersect(names(control$outOptimControl), optimControlArgNames)
      numValidNames <- length(validNames)
      if(numValidNames > 0){
        for(i in 1:numValidNames){
          outOptControl[[validNames[i]]] <- control$outOptimControl[[validNames[i]]]
        }   
      }
    }
    if(!is.null(control$innerOptimControl)) {
      validNames_inner <- intersect(names(control$innerOptimControl), optimControlArgNames)
      numValidNames_inner <- length(validNames_inner)
      if(numValidNames_inner > 0){
        for(i in 1:numValidNames_inner)
          innerOptControl[[validNames_inner[i]]] <- control$innerOptimControl[[validNames_inner[i]]]   
      }
    }
    outOptControl$fnscale <- -1 
    innerOptControl$fnscale <- -1 
    if(!is.null(control$innerOptimMethod) && (control$innerOptimMethod %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B"))){
      innerOptMethod <- control$innerOptimMethod
    }
    else innerOptMethod <- "BFGS"
    
    ## Create an AGHQuad (Adaptive Gauss-Hermite Quadrature) nimbleFunctionList
    AGHQuad_nfl <- nimbleFunctionList(AGHQuad_BASE)
    scalarRENodes <- model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE)
    nre <- length(scalarRENodes)
    if(nre > 0){
      ## Record the order of random effects processed internally
      internalRandomEffectsNodes <- NULL
      lenInternalRENodeSets <- NULL
      if(isFALSE(split)) { ## Do all randomEffectsNodes in one set
        internalRandomEffectsNodes <- randomEffectsNodes
        lenInternalRENodeSets <- nre
        if(is.null(control$innerOptimStart)) innerOptStart <- values(model, randomEffectsNodes)
        else {
          providedStart <- control$innerOptimStart
          if(any(providedStart %in% c("last", "last.best"))) innerOptStart <- providedStart
          else if(is.numeric(sum(providedStart)) && (length(providedStart) == nre)) innerOptStart <- providedStart
          else innerOptStart <- values(model, randomEffectsNodes)
        }
        ## In case random effects are not properly initialized
        if(!any(innerOptStart %in% c("last", "last.best")) & any(is.infinite(innerOptStart) | is.na(innerOptStart) | is.nan(innerOptStart))){
          all_reTransform <- parameterTransform(model, randomEffectsNodes)
          all_reTransform_length <- all_reTransform$getTransformedLength()
          innerOptStart <- all_reTransform$inverseTransform(rep(0, all_reTransform_length))
        }
        ## Build AGHQuad
        if(nre > 1) AGHQuad_nfl[[1]] <- buildOneAGHQuad(model, paramNodes, randomEffectsNodes, calcNodes, innerOptControl, innerOptMethod, innerOptStart)
        else AGHQuad_nfl[[1]] <- buildOneAGHQuad1D(model, paramNodes, randomEffectsNodes, calcNodes, innerOptControl, "CG", innerOptStart)
      }
      else {## Split randomEffectsNodes into conditionally independent sets
        reSets <- MargNodes$randomEffectsSets
        num_reSets <- length(reSets)
        if(num_reSets == 0){
          stop("There was a problem determining conditionally independent random effects sets for this model.")
        }
        for(i in seq_along(reSets)){
          ## Work with one conditionally independent set of latent states
          these_reNodes <- reSets[[i]]
          internalRandomEffectsNodes <- c(internalRandomEffectsNodes, these_reNodes)
          ## find paramNodes and calcNodes for this set of reNodes
          ## paramNodes are the same for all AGHQuad_nfl elements. In the future this could be customized.
          these_reDeps <- model$getDependencies(these_reNodes)  ## candidate calcNodes via reNodes
          these_calcNodes <- intersect(calcNodes, these_reDeps) ## definite calcNodes
          nre_these <- length(model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE))
          lenInternalRENodeSets <- c(lenInternalRENodeSets, nre_these)
          ## Process start values for inner optimisation
          if(is.null(control$innerOptimStart)) innerOptStart <- values(model, these_reNodes)
          else {
            providedStart <- control$innerOptimStart
            if(any(providedStart %in% c("last", "last.best"))) innerOptStart <- providedStart
            else if(is.numeric(sum(providedStart)) && (length(providedStart) == nre)){
              these_reNodes_inds <- unlist(lapply(model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE), function(x) {which(scalarRENodes == x)}))
              innerOptStart <- providedStart[these_reNodes_inds]
            }
            else innerOptStart <- values(model, these_reNodes)
          }
          ## In case random effects are not properly initialized
          if(!any(innerOptStart %in% c("last", "last.best")) & any(is.infinite(innerOptStart) | is.na(innerOptStart) | is.nan(innerOptStart))){
            these_reTransform <- parameterTransform(model, these_reNodes)
            these_reTransform_length <- these_reTransform$getTransformedLength()
            innerOptStart <- these_reTransform$inverseTransform(rep(0, these_reTransform_length))
          }
          ## Build AGHQuad for each set
          if(nre_these > 1){
            AGHQuad_nfl[[i]] <- buildOneAGHQuad(model, paramNodes, these_reNodes, these_calcNodes, innerOptControl, innerOptMethod, innerOptStart)
          }
          else AGHQuad_nfl[[i]] <- buildOneAGHQuad1D(model, paramNodes, these_reNodes, these_calcNodes, innerOptControl, "CG", innerOptStart)
        }
      }
      if(length(lenInternalRENodeSets) == 1) lenInternalRENodeSets <- c(lenInternalRENodeSets, -1)
      reTransform <- parameterTransform(model, internalRandomEffectsNodes)
      reTransform_length <- reTransform$getTransformedLength()
      if(reTransform_length > 1) reTransform_indices <- 1:reTransform_length
      else reTransform_indices <- c(1, -1)
      
      reNodesAsScalars <- model$expandNodeNames(internalRandomEffectsNodes, returnScalarComponents = TRUE)
      reNodesAsScalars_vec <- reNodesAsScalars
      if(nre == 1) reNodesAsScalars_vec <- c(reNodesAsScalars, "_EXTRA_")
      reNodesAsScalars_first <- reNodesAsScalars[1]
    }
    else{
      ## No random effects
      lenInternalRENodeSets <- numeric(2)
      reTransform <- parameterTransform(model, paramNodes[1], control = list(allowDeterm = allowNonPriors)) ## Won't be needed at all
      reTransform_indices <- numeric(2)
      reNodesAsScalars_vec <- character(0)
      reNodesAsScalars_first <- character(1)
      if(num_calcNodesOther == 0)
        stop("Both calcNodesOther and randomEffectsNodes are empty for Laplace or AGHQuad for the given model.")
    }
    
    paramNodesAsScalars     <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
    npar <- length(paramNodesAsScalars)
    paramNodesAsScalars_vec <- paramNodesAsScalars
    if(npar == 1) paramNodesAsScalars_vec <- c(paramNodesAsScalars, "_EXTRA_")
    paramNodesAsScalars_first <- paramNodesAsScalars[1]
    if(npar == 1) p_indices <- c(1, -1)
    else p_indices <- 1:npar
    ## setupOutputs(reNodesAsScalars, paramNodesAsScalars)
    
    ## Automated transformation for parameters
    paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = allowNonPriors))
    pTransform_length <- paramsTransform$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
    
    ## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    ## Default calculation method for AGHQuad
    methodID <- 2
    ## The nimbleList definitions AGHQuad_params and AGHQuad_summary
    ## have moved to predefined nimbleLists.
  },## End of setup
  run = function(){},
  methods = list(
    get_node_names_vec = function(returnParams = logical(0, default = TRUE)) {
      returnType(character(1))
      if(returnParams) return(paramNodesAsScalars_vec)
      else return(reNodesAsScalars_vec)
    },
    get_node_names_single = function(returnParams = logical(0, default = TRUE)) {
      returnType(character())
      if(returnParams) return(paramNodesAsScalars_first)
      else return(reNodesAsScalars_first)
    },
    set_method = function(method = integer()) {
      if(nre == 0) print("AGHQuad or Laplace approximation is not needed for the given model: no random effects")
      if(!any(c(1, 2, 3) == method)) stop("Choose a valid method ID from 1, 2, and 3")
      methodID <<- method
    },
    get_method = function() {
      return(methodID)
      returnType(integer())
    },
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(pTransform_length == 1){
        if(length(pTransform_indices) == 2){
          pTransform_indices <<- numeric(length = 1, value = 1)
        }
      }
      if(npar == 1){
        if(length(p_indices) == 2){
          p_indices <<- numeric(length = 1, value = 1)
        }
      }
      one_time_fixes_done <<- TRUE
    },
    ## Other log-likelihood (parts not involving random effects, i.e. simply
    ## additional calculations in the model) in terms of original parameters
    otherLogLik = function(p = double(1)) {
      if(num_calcNodesOther == 0) stop("calcNodesOther is empty: there is no exact likelihood component for the model")
      values(model, paramNodes) <<- p
      ans <- model$calculate(calcNodesOther)
      return(ans)
      returnType(double())
    },
    ## Gradient of the exact log-likelihood w.r.t parameters
    gr_otherLogLik_internal = function(p = double(1)) {
      if(num_calcNodesOther == 0) stop("calcNodesOther is empty: there is no exact likelihood component for the model")
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(otherLogLik(p), wrt = p_indices, order = 1, model = model,
                    updateNodes = otherLogLik_updateNodes, constantNodes = otherLogLik_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    gr_otherLogLik = function(p = double(1)) {
      if(num_calcNodesOther == 0) stop("calcNodesOther is empty: there is no exact likelihood component for the model")
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(gr_otherLogLik_internal(p), wrt = p_indices, order = 0, model = model,
                    updateNodes = otherLogLik_updateNodes, constantNodes = otherLogLik_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## AGHQuad approximation in terms of original parameters
    calcLogLik = function(p = double(1), trans = logical(0, default = FALSE)) {
      if(!one_time_fixes_done) one_time_fixes()
      if(trans) {
        if(length(p) != pTransform_length) {
          print("  [Warning] For calcLogLik (or calcLaplace) with trans = TRUE, p should be length ", pTransform_length, " but was provided with length ", length(p),".")
          stop("Wrong length for p in calcLogLik  (or calcLaplace) with trans = TRUE.")
        }
        p <- paramsTransform$inverseTransform(p)
      }
      if(length(p) != npar) {
        print("  [Warning] For calcLogLik (or calcLaplace), p should be length ", npar, " but is length ", length(p), ".")
        stop("Wrong length for p in calcLogLik  (or calcLaplace).")
      }
      if(num_calcNodesOther > 0) ans <- otherLogLik(p)
      else ans <- 0
      if(nre > 0){
        for(i in seq_along(AGHQuad_nfl)){
          if(methodID == 1) ans <- ans + AGHQuad_nfl[[i]]$calcLogLik1(p)
          else if(methodID == 2) ans <- ans + AGHQuad_nfl[[i]]$calcLogLik2(p)
          else ans <- ans + AGHQuad_nfl[[i]]$calcLogLik3(p)
        }
      }
      if(is.nan(ans) | is.na(ans)) ans <- -Inf
      return(ans)
      returnType(double())
    },
    calcLaplace = function(p = double(1), trans = logical(0, default = FALSE)) {
      ans <- calcLogLik(p, trans)
      return(ans)
      returnType(double())
    },
    ## Gradient of the AGHQuad approximation w.r.t. parameters
    gr_logLik = function(p = double(1), trans = logical(0, default=FALSE)) {
      if(!one_time_fixes_done) one_time_fixes()
      if(trans) {
        if(length(p) != pTransform_length) {
          print("  [Warning] For gr_logLik (or gr_Laplace) with trans = TRUE, p should be length ", pTransform_length, " but was provided with length ", length(p),".")
          stop("Wrong length for p in gr_logLik (or gr_Laplace) with trans = TRUE.")
        }
        pDerivs <- derivs_pInverseTransform(p, c(0, 1))
        p <- pDerivs$value
      }
      if(length(p) != npar) {
        print("    [Warning] For gr_logLik (or gr_Laplace), p should be length ", npar, " but is length ", length(p), ".")
        stop("Wrong length for p in gr_logLik (or gr_Laplace).")
      }
      if(num_calcNodesOther > 0) ans <- gr_otherLogLik(p)
      else ans <- numeric(length = npar)
      if(nre > 0){
        for(i in seq_along(AGHQuad_nfl)){
          if(methodID == 1) ans <- ans + AGHQuad_nfl[[i]]$gr_logLik1(p)
          else if(methodID == 2) ans <- ans + AGHQuad_nfl[[i]]$gr_logLik2(p)
          else ans <- ans + AGHQuad_nfl[[i]]$gr_logLik3(p)
        }
      }
      if(trans) {
        ans <- (ans %*% pDerivs$jacobian)[1,]
      }
      return(ans)
      returnType(double(1))
    },
    gr_Laplace = function(p = double(1), trans = logical(0, default=FALSE)) {
      ans <- gr_logLik(p, trans)
      return(ans)
      returnType(double(1))
    },
    ## AGHQuad approximation in terms of transformed parameters
    calcLogLik_pTransformed = function(pTransform = double(1)) {
      ans <- calcLogLik(pTransform, trans = TRUE)
      ## if(!one_time_fixes_done) one_time_fixes()
      ## p <- paramsTransform$inverseTransform(pTransform)
      ## ans <- calcLogLik(p)
      ## if(is.nan(ans) | is.na(ans)) ans <- -Inf
      return(ans)
      returnType(double())
    },
    ## Inverse transform parameters to original scale
    pInverseTransform = function(pTransform = double(1)) {
      p <- paramsTransform$inverseTransform(pTransform)
      return(p)
      returnType(double(1))
    },
    ## Jacobian of the inverse transformation for parameters
    derivs_pInverseTransform = function(pTransform = double(1), order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(pInverseTransform(pTransform), wrt = pTransform_indices, order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    ## Inverse transform random effects to original scale
    reInverseTransform = function(reTrans = double(1)) {
      if(nre == 0) stop("No random effects in the model")
      re <- reTransform$inverseTransform(reTrans)
      return(re)
      returnType(double(1))
    },
    ## Jacobian of the inverse transformation
    derivs_reInverseTransform = function(reTrans = double(1), order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      if(nre == 0) stop("No random effects in the model")
      ans <- derivs(reInverseTransform(reTrans), wrt = reTransform_indices, order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    ## Gradient of the AGHQuad approximation in terms of transformed parameters
    gr_logLik_pTransformed = function(pTransform = double(1)) {
      ans <- gr_logLik(pTransform, trans = TRUE)
      ## if(!one_time_fixes_done) one_time_fixes()
      ## pDerivs <- derivs_pInverseTransform(pTransform, c(0, 1))
      ## gr <- gr_logLik(pDerivs$value) ## pDerivs$value gives original param values
      ## ans <- (gr %*% pDerivs$jacobian)[1,]
      return(ans)
      returnType(double(1))
    },
    ## Calculate MLEs of parameters
    findMLE = function(pStart  = double(1, default = Inf),
                       method  = character(0, default = "BFGS"),
                       hessian = logical(0, default = TRUE)) {
      if(any(abs(pStart) == Inf)) pStart <- values(model, paramNodes)
      if(length(pStart) != npar) {
        print("  [Warning] For findMLE, pStart should be length ", npar, " but is length ", length(pStart), ".")
        ans <- optimResultNimbleList$new()
        return(ans)
#        stop("Wrong length for pStart in findMLE.")
      }
      ## In case parameter nodes are not properly initialized
      if(any_na(pStart) | any_nan(pStart) | any(abs(pStart)==Inf)) pStartTransform <- rep(0, pTransform_length)
      else pStartTransform <- paramsTransform$transform(pStart)
      ## In case bad start values are provided 
      if(any_na(pStartTransform) | any_nan(pStartTransform) | any(abs(pStartTransform)==Inf)) pStartTransform <- rep(0, pTransform_length)
      optRes <- optim(pStartTransform, calcLogLik_pTransformed, gr_logLik_pTransformed, method = method, control = outOptControl, hessian = hessian)
      if(optRes$convergence != 0) 
        print("Warning: optim has a non-zero convergence code: ", optRes$convergence, ".\n",
              "The control parameters of optim can be adjusted in the control argument of\n",
              "buildLaplace or buildAGHQuad via list(outOptimControl = list()).")
      ## Back transform results to original scale
      optRes$par <- paramsTransform$inverseTransform(optRes$par)
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Optimized random effects given transformed parameter values
    optimRandomEffects = function(pTransform = double(1)){
      if(nre == 0) stop("No random effects in the model")
      p <- pInverseTransform(pTransform)
      raneff <- numeric(nre)
      tmp <- numeric(nre) ## Not sure this is needed. 
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        if(methodID == 1) tmp <- AGHQuad_nfl[[i]]$update_max_inner_logLik_internal(p)
        else tmp <- AGHQuad_nfl[[i]]$update_max_inner_logLik(p)
        numre <- dim(tmp)[1]
        raneff[(tot+1):(tot+numre)] <- tmp
        tot <- tot + numre
      }
      return(raneff)
      returnType(double(1))
    },
    ## Inverse of the negative Hessian of log-likelihood wrt transformed random effects
    inverse_negHess = function(p = double(1), reTransform = double(1)){
      if(nre == 0) stop("No random effects in the model")
      invHess <- matrix(value = 0, nrow = nre, ncol = nre)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        numre <- lenInternalRENodeSets[i]
        tmp <- AGHQuad_nfl[[i]]$negHess(p, reTransform[(tot+1):(tot+numre)])
        invHess[(tot+1):(tot+numre), (tot+1):(tot+numre)] <- inverse(tmp)
        tot <- tot + numre
      }
      return(invHess)
      returnType(double(2))
    },
    ## Hessian of joint log-likelihood wrt parameters and (transformed) random effects
    hess_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)){
      if(nre == 0) stop("No random effects in the model")
      ans <- matrix(value = 0, nrow = npar, ncol = nre)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        numre <- lenInternalRENodeSets[i]
        if(methodID == 1) tmp <- AGHQuad_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform[(tot+1):(tot+numre)])
        else tmp <- AGHQuad_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re(p, reTransform[(tot+1):(tot+numre)])
        ans[1:npar, (tot+1):(tot+numre)] <- tmp
        tot <- tot + numre
      }
      return(ans)
      returnType(double(2))
    },
    ## Summarise AGHQuad MLE results
    summary = function(MLEoutput          = optimResultNimbleList(),
                       originalScale             = logical(0, default = TRUE),
                       calcRandomEffectsStdError = logical(0, default = FALSE), 
                       returnJointCovariance     = logical(0, default = FALSE)){
      if(dim(MLEoutput$hessian)[1] == 0) stop("Hessian matrix was not calculated for Laplace or AGHQuad MLE")
      ## Output lists
      ans <- AGHQuad_summary$new()
      pres <- AGHQuad_params$new()
      ranres <- AGHQuad_params$new()
      ## Parameters
      p <- MLEoutput$par
      pTransform <- paramsTransform$transform(p)
      vcov_pTransform <- -inverse(MLEoutput$hessian)
      stdErr_pTransform <- sqrt(diag(vcov_pTransform))
      if(nre == 0) { ## No random effects
        ranres$estimates <- numeric(0)
        ranres$stdErrors <- numeric(0)
        if(originalScale){
          derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
          JacobpInvTransform   <- derivspInvTransform$jacobian
          stdErr_p <- numeric(npar)
          if(returnJointCovariance) {
            vcov <- JacobpInvTransform %*% vcov_pTransform %*% t(JacobpInvTransform)
            stdErr_p <- sqrt(diag(vcov))
            ans$vcov <- vcov
          }
          else{
            for(i in 1:npar){
              var_p_i <- (JacobpInvTransform[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobpInvTransform[i,,drop=FALSE]))[1,1]
              stdErr_p[i] <- sqrt(var_p_i)
            }
            ans$vcov <- matrix(nrow = 0, ncol = 0)
          }
          pres$estimates <- p
          pres$stdErrors <- stdErr_p
        }
        else {
          pres$estimates <- pTransform
          pres$stdErrors <- stdErr_pTransform
          if(returnJointCovariance) ans$vcov <- vcov_pTransform
          else ans$vcov <- matrix(0, nrow = 0, ncol = 0)
        }
      }
      else{
        ## Random effects
        optreTransform <- optimRandomEffects(pTransform)
        optre <- reInverseTransform(optreTransform)
        ntot <- npar + nre
        if(returnJointCovariance) {
          ## Inverse of the negative Hessian of log-likelihood wrt transformed random effects at MLEs
          inv_negHess <- inverse_negHess(p, optreTransform)
          jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
          jointInvNegHessZero[1:nre, 1:nre] <- inv_negHess
          ## Hessian of log-likelihood wrt to params and transformed random effects
          hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
          ## Derivative of inverse transformation for params
          derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
          JacobpInvTransform   <- derivspInvTransform$jacobian
          ## Jacobian of optimized random effects wrt transformed parameters
          JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
          jointJacob <- matrix(NA, nrow = ntot, ncol = npar)
          jointJacob[1:nre, 1:npar] <- JacobOptreWrtParams
          jointJacob[(nre+1):ntot, 1:npar] <- diag(npar)
          ## Joint covariance matrix on transformed scale
          vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
          if(originalScale){
            derivs_reInvTransform <- derivs_reInverseTransform(optreTransform, c(0, 1))
            Jacob_reInvTransform  <- derivs_reInvTransform$jacobian
            Jacob_JointInvTransform <- matrix(0, nrow = ntot, ncol = ntot)
            Jacob_JointInvTransform[1:nre, 1:nre] <- Jacob_reInvTransform
            Jacob_JointInvTransform[(nre+1):ntot, (nre+1):ntot] <- JacobpInvTransform
            vcov <- Jacob_JointInvTransform %*% vcov_Transform %*% t(Jacob_JointInvTransform)
            stdErr_re_p <- sqrt(diag(vcov))
            stdErr_p <- stdErr_re_p[(nre+1):ntot]
            if(calcRandomEffectsStdError){
              ranres$stdErrors <- stdErr_re_p[1:nre]
            }
            else{
              ranres$stdErrors <- numeric(0)
            }
            ans$vcov <- vcov
            pres$estimates <- p
            pres$stdErrors <- stdErr_p
            ranres$estimates <- optre
          }## End of if(originalScale)
          else { ## On transformed scale
            if(calcRandomEffectsStdError){
              stdErr_reTransform <- sqrt(diag(vcov_Transform)[1:nre])
              ranres$stdErrors <- stdErr_reTransform
            }
            else{
              ranres$stdErrors <- numeric(0)
            }
            ans$vcov <- vcov_Transform
            pres$estimates <- pTransform
            pres$stdErrors <- sqrt(diag(vcov_Transform)[(nre+1):ntot])
            ranres$estimates <- optreTransform
          }
        }## End of if(returnJointCovariance)
        else { ## Do not return joint covariance matrix
          ans$vcov <- matrix(nrow = 0, ncol = 0)
          if(originalScale){## On original scale
            pres$estimates <- p
            ranres$estimates <- optre
            if(calcRandomEffectsStdError){
              ## Joint covariance matrix on transform scale
              inv_negHess <- inverse_negHess(p, optreTransform)
              jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
              jointInvNegHessZero[1:nre, 1:nre] <- inv_negHess
              ## Hessian of log-likelihood wrt to params and transformed random effects
              hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
              ## Derivative of inverse transformation for params
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              ## Jacobian of optimized random effects wrt transformed parameters
              JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
              jointJacob <- matrix(NA, nrow = ntot, ncol = npar)
              jointJacob[1:nre, 1:npar] <- JacobOptreWrtParams
              jointJacob[(nre+1):ntot, 1:npar] <- diag(npar)
              ## Join covariance matrix on transformed scale
              vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
              ## Derivatives information
              derivs_reInvTransform <- derivs_reInverseTransform(optreTransform, c(0, 1))
              Jacob_reInvTransform  <- derivs_reInvTransform$jacobian
              Jacob_JointInvTransform <- matrix(0, nrow = ntot, ncol = ntot)
              Jacob_JointInvTransform[1:nre, 1:nre] <- Jacob_reInvTransform
              Jacob_JointInvTransform[(nre+1):ntot, (nre+1):ntot] <- JacobpInvTransform
              stdErr <- numeric(ntot)
              for(i in 1:ntot){
                var_i <- (Jacob_JointInvTransform[i,,drop=FALSE] %*% vcov_Transform %*% t(Jacob_JointInvTransform[i,,drop=FALSE]))[1,1]
                stdErr[i] <- sqrt(var_i)
              }
              stdErr_p <- stdErr[(nre+1):ntot]
              stdErr_re <- stdErr[1:nre]
              pres$stdErrors   <- stdErr_p
              ranres$stdErrors <- stdErr_re
            }## End of if(calcRandomEffectsStdError)
            else { ## Do not calculate standard errors of random effects estimates
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              stdErr_p <- numeric(npar)
              for(i in 1:npar){
                var_p_i <- (JacobpInvTransform[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobpInvTransform[i,,drop=FALSE]))[1,1]
                stdErr_p[i] <- sqrt(var_p_i)
              }
              pres$stdErrors <- stdErr_p
              ranres$stdErrors <- numeric(0)
            }
          }## End of if(originalScale)
          else {## On transformed scale
            pres$estimates <- pTransform
            pres$stdErrors <- stdErr_pTransform
            ranres$estimates <- optreTransform
            if(calcRandomEffectsStdError){
              inv_negHess <- inverse_negHess(p, optreTransform)
              jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
              jointInvNegHessZero[1:nre, 1:nre] <- inv_negHess
              ## Hessian of log-likelihood wrt to params and transformed random effects
              hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
              ## Derivative of inverse transformation for params
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              ## Jacobian of optimized random effects wrt transformed parameters
              JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
              stdErr_reTransform <- numeric(nre)
              for(i in 1:nre){
                var_reTransform_i <- inv_negHess[i, i] + (JacobOptreWrtParams[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobOptreWrtParams[i,,drop=FALSE]))[1,1]
                stdErr_reTransform[i] <- sqrt(var_reTransform_i)
              }
              ranres$stdErrors <- stdErr_reTransform
            }
            else{
              ranres$stdErrors <- numeric(0)
            }
          }
        }
      }
      pres$names <- paramNodesAsScalars_vec
      ranres$names <- reNodesAsScalars_vec
      ans$params <- pres
      ans$randomEffects <- ranres
      if(originalScale) ans$scale <- "original"
      else ans$scale <- "transformed"
      return(ans)
      returnType(AGHQuad_summary())
    }
  ),
  buildDerivs = list(pInverseTransform  = list(),
                     reInverseTransform = list(),
                     otherLogLik = list(),
                     gr_otherLogLik_internal = list()
                     )
)

#' Laplace approximation
#' 
#' Builds a Laplace approximation algorithm for a given NIMBLE model. 
#' 
#' @param model a NIMBLE model object, such as returned by \code{nimbleModel}.
#'   The model must have automatic derivatives (AD) turned on, e.g. by using
#'   \code{buildDerivs=TRUE} in \code{nimbleModel}.
#' @param paramNodes a character vector of names of parameter nodes in the
#'   model; defaults to top-level stochastic nodes, as determined by
#'   \code{model$getNodeNames(topOnly=TRUE, stochOnly=TRUE)}. Alternatively,
#'   \code{paramNodes} can be a list in the format returned by
#'   \code{setupMargNodes}, in which case \code{randomEffectsNodes},
#'   \code{calcNodes}, and \code{calcNodesOther} are not needed (and will be
#'   ignored).
#' @param randomEffectsNodes a character vector of names of unobserved (latent)
#'   nodes to marginalize (integrate) over using Laplace approximation; defaults
#'   to latent stochastic nodes that depend on \code{paramNodes}.
#' @param calcNodes a character vector of names of nodes for calculating the
#'   integrand for Laplace approximation; defaults to nodes that depend on
#'   \code{randomEffectsNodes}, as determined by
#'   \code{model$geteDependencies(randomEffectsNodes)} (which will include
#'   \code{randomEffectsNodes}). There may be deterministic nodes between
#'   \code{paramNodes} and \code{randomEffectsNodes}. These will be included in
#'   calculations automatically and thus do not need to be included in
#'   \code{calcNodes} (but there is no problem if they are).
#' @param calcNodesOther a character vector of names of nodes for calculating
#'   terms in the log-likelihood that do not depend on any
#'   \code{randomEffectsNodes}, and thus are not part of the marginalization,
#'   but should be included for purposes of finding the MLE. This defaults to
#'   stochastic nodes that depend on \code{paramNodes} but are not part of and
#'   do not depend on \code{randomEffectsNodes}. There may be deterministic
#'   nodes between \code{paramNodes} and \code{calcNodesOther}. These will be
#'   included in calculations automatically and thus do not need to be included
#'   in \code{calcNodesOther} (but there is no problem if they are).
#' @param control a named list for providing
#'   additional settings used in Laplace approximation. See
#'   \code{control} section below.
#'
#' @section \code{buildLaplace}:
#'
#' \code{buildLaplace} is the main function for constructing the Laplace
#' approximation for a given model.
#'
#' One only needs to provide a NIMBLE model object and then the function will
#' construct the pieces necessary for Laplace approximation to marginalize over
#' all latent states (aka random effects) in a model. To do so, it will
#' determine default values for \code{paramNodes}, \code{randomEffectsNodes},
#' \code{calcNodes}, and \code{calcNodesOther} as described above.
#'
#' The default values are obtained by calling \code{setupMargNodes}, whose
#' arguments match those here (except for a few arguments which are taken from
#' control list elements here). One can call that function to see exactly how
#' nodes will be arranged for Laplace approximation. One can also call it,
#' customize the returned list, and then provide that to \code{buildLaplace} as
#' \code{paramNodes}.
#'
#' If any \code{paramNodes} (parameters) or \code{randomEffectsNodes} (random
#' effects / latent states) have constraints on the range of valid values
#' (because of the distribution they follow), they will be used on a transformed
#' scale determined by \code{parameterTransform}. This means the Laplace
#' approximation itself will be done on the transformed scale for random effects
#' and finding the MLE will be done on the transformed scale for parameters. For
#' parameters, any prior distributions are not included in calculations, but
#' they are used to determine valid parameter ranges. For example, if
#' \code{sigma} is a standard deviation, declare it with a prior such as
#' \code{sigma ~ dhalfflat()} to indicate that it must be greater than 0.
#'
#' For default determination of parameters, all parameters must have a prior
#' distribution simply to indicate the range of valid values. For a param
#' \code{p} that has no constraint, a simple choice is \code{p ~ dflat()}.
#'
#' The object returned by \code{buildLaplace} is a nimbleFunction object with
#' numerous methods (functions). The most useful ones are:
#'
#' \itemize{
#'
#' \item \code{calcLogLik(p, trans)}. Calculate the Laplace approximation to
#'       the marginal log-likelihood function at parameter value \code{p}, which
#'       (if \code{trans} is FALSE, which is the default) should match the order
#'       of \code{paramNodes}. For any non-scalar nodes in \code{paramNodes},
#'       the order within the node is column-major (which can be seen for R
#'       objects using \code{as.numeric}). Return value is the scalar
#'       (approximate, marginal) log likelihood.
#'
#'       If \code{trans} is TRUE, then \code{p} is the vector of parameters on
#'       the transformed scale, if any, described above. In this case, the
#'       parameters on the original scale (as the model was written) will be
#'       determined by calling the method \code{pInverseTransform(p)}. Note that
#'       the length of the parameter vector on the transformed scale might not
#'       be the same as on the original scale (because some constraints of
#'       non-scalar parameters result in fewer free transformed parameters than
#'       original parameters).
#'
#' \item \code{calcLaplace(p, trans)}. This is the same as \code{calcLaplace}.
#'
#' \item \code{findMLE(pStart, method, hessian)}. Find the maximum likelihood
#'         estimates of the Laplace-approximated marginal likelihood. Arguments
#'         include \code{pStart}: initial parameter values (defaults to
#'         parameter values currently in the model); \code{method}: (outer)
#'         optimization method to use in \code{optim} (defaults to "BFGS"); and
#'         \code{hessian}: whether to calculate and return the Hessian matrix
#'         (defaults to \code{TRUE}). Second derivatives in the Hessian are
#'         determined by finite differences of the gradients obtained by
#'         automatic differentiation (AD). Return value is a nimbleList of type
#'         \code{optimResultNimbleList}, similar to what is returned by R's
#'         optim.
#'
#' \item \code{summary(MLEoutput, originalScale,
#'        calcRandomEffectsStdError, returnJointCovariance)}. Summarize the
#'        maximum likelihood estimation results, given object
#'        \code{MLEoutput} that was returned by \code{findMLE}. The
#'        summary can include a covariance matrix for the parameters, the random
#'        effects, or both, and these can be returned on the original parameter
#'        scale or on the (potentially) transformed scale(s) used in estimation.
#'
#' In addition, \code{summary} accepts the following optional arguments:
#'
#'        \itemize{
#'
#'           \item \code{originalScale}. Logical. If TRUE, the function returns
#'           results on the original scale(s) of parameters and random effects;
#'           otherwise, it returns results on the transformed scale(s). If there
#'           are no constraints, the two scales are identical. Defaults to TRUE.
#'
#'           \item \code{calcRandomEffectsStdError}. Logical. If TRUE, standard
#'           errors of random effects will be calculated.
#'           Defaults to FALSE.
#'
#'           \item \code{returnJointCovariance}. Logical. If TRUE, the joint
#'           variance-covariance matrix of the parameters and the random effects
#'           will be returned. Defaults to FALSE.
#'
#'        }
#'
#'        The object returned by \code{summary} is a nimbleList with elements:
#'
#'        \itemize{
#'
#'           \item \code{params}. A list that contains estimates and standard
#'           errors of parameters (on the original or transformed scale, as
#'           chosen by \code{originalScale}).
#'
#'           \item \code{random}. A list that contains estimates of random
#'           effects and, if requested (\code{calcRandomEffectsStdError=TRUE})
#'           their standard errors, on original or transformed scale. Standard
#'           errors are calculated following the generalized delta method of
#'           Kass and Steffey (1989).
#'
#'           \item \code{vcov}. If requested (i.e.
#'           \code{returnJointCovariance=TRUE}), the joint variance-covariance
#'           matrix of the random effects and parameters, on original or
#'           transformed scale.
#'
#'           \item \code{scale}. \code{"original"} or \code{"transformed"}, the
#'        scale on which results were requested.
#'
#'        }
#'
#'     }
#'
#' Additional methods to access or control more details of the Laplace approximation include:
#'
#' \itemize{
#'
#'   \item \code{get_node_name_vec(returnParams)}. Return a vector (>1) of names
#'   of parameters/random effects nodes, according to \code{returnParams =
#'   TRUE/FALSE}. Use this if there is more than one parameter.
#'
#'   \item \code{get_node_name_single(returnParams)}. Return the name of a
#'   single parameter/random effect node, \code{returnParams = TRUE/FALSE}. Use
#'   this if there is only one parameter.
#'
#'   \item \code{set_method(method)}. Set method ID for calculating the Laplace
#'   approximation and gradient: 1 (\code{Laplace1}), 2 (\code{Laplace2},
#'   default method), or 3 (\code{Laplace3}). See below for more details. Users
#'   wanting to explore efficiency can try switching from method 2 (default) to
#'   methods 1 or 3 and comparing performance. The first Laplace approximation
#'   with each method will be (much) slower than subsequent Laplace
#'   approximations.
#'
#'   \item \code{get_method()}. Return the current method ID.
#'
#'   \item \code{gr_logLik(p, trans)}. Gradient of the Laplace-approximated
#'   marginal likelihood at parameter value \code{p}. Argument \code{trans} is
#'   similar to that in \code{calcLaplace}.
#'
#'   \item \code{gr_Laplace(p, trans)}. This is the same as \code{gr_logLik}.
#'
#'   \item \code{otherLogLik(p)}. Calculate the \code{calcNodesOther}
#'   nodes, which returns the log-likelihood of the parts of the model that are
#'   not included in the Laplace approximation. \code{p} is the vector of
#'   parameter values in the order of \code{paramNames}.
#'
#'   \item \code{gr_otherLogLik(p)}. Gradient (vector of derivatives with
#'   respect to each parameter) of \code{otherLogLik(p)}. This is obtained
#'   using automatic differentiation (AD) with double-taping. Results should
#'   match \code{gr_otherLogLik_internal(p)} but may be more efficient after
#'   the first call.
#'
#' }
#'
#' Finally, methods that are primarily for internal use by other methods include:
#'
#' \itemize{
#'
#'    \item \code{p_transformed_gr_Laplace(pTransform)}. Gradient of the Laplace
#'     approximation (\code{p_transformed_Laplace(pTransform)}) at transformed 
#'     (unconstrained) parameter value \code{pTransform}.
#'
#'    \item \code{pInverseTransform(pTransform)}. Back-transform the transformed
#'    parameter value \code{pTransform} to original scale.
#'
#'    \item \code{derivs_pInverseTransform(pTransform, order)}. Derivatives of
#'    the back-transformation (i.e. inverse of parameter transformation) with
#'    respect to transformed parameters at \code{pTransform}. Derivative order
#'    is given by \code{order} (any of 0, 1, and/or 2).
#'
#'    \item \code{reInverseTransform(reTrans)}. Back-transform the transformed
#'    random effects value \code{reTrans} to original scale.
#'
#'    \item \code{derivs_reInverseTransform(reTrans, order)}. Derivatives of the
#'    back-transformation (i.e. inverse of random effects transformation) with
#'    respect to transformed random effects at \code{reTrans}. Derivative order
#'    is given by \code{order} (any of 0, 1, and/or 2).
#'
#'    \item \code{optimRandomEffects(pTransform)}. Calculate the optimized
#'    random effects given transformed parameter value \code{pTransform}. The
#'    optimized random effects are the mode of the conditional distribution of
#'    random effects given data at parameters \code{pTransform}, i.e. the
#'    calculation of \code{calcNodes}.
#'
#'    \item \code{inverse_negHess(p, reTransform)}. Calculate the inverse of the
#'    negative Hessian matrix of the joint (parameters and random effects)
#'    log-likelihood with respect to transformed random effects, evaluated at
#'    parameter value \code{p} and transformed random effects
#'    \code{reTransform}.
#'
#'    \item \code{hess_logLik_wrt_p_wrt_re(p, reTransform)}. Calculate the
#'    Hessian matrix of the joint log-likelihood with respect to parameters and
#'    transformed random effects, evaluated at parameter value \code{p} and
#'    transformed random effects \code{reTransform}.
#'
#'   \item \code{one_time_fixes()}. Users do not need to run this. Is is called
#'   when necessary internally to fix dimensionality issues if there is only
#'   one parameter in the model.
#'
#'   \item \code{p_transformed_Laplace(pTransform)}. Laplace approximation at
#'         transformed (unconstrained) parameter value \code{pTransform}. To
#'         make maximizing the Laplace likelihood unconstrained, an automated
#'         transformation via \code{\link{parameterTransform}} is performed on
#'         any parameters with constraints indicated by their priors (even
#'         though the prior probabilities are not used).
#'
#'   \item \code{gr_otherLogLik_internal(p)}. Gradient (vector of
#'   derivatives with respect to each parameter) of \code{otherLogLik(p)}.
#'   This is obtained using automatic differentiation (AD) with single-taping.
#'   First call will always be slower than later calls.
#'
#' }
#'
#' @section \code{control} list:
#' 
#' \code{buildLaplace} accepts the following control list elements:
#'
#' \itemize{
#'
#'   \item \code{split}. If TRUE (default), \code{randomEffectsNodes} will be
#'         split into conditionally independent sets if possible. This
#'         facilitates more efficient Laplace approximation because each
#'         conditionally independent set can be marginalized independently. If
#'         FALSE, \code{randomEffectsNodes} will be handled as one multivariate
#'         block, with one multivariate Laplace approximation. If \code{split}
#'         is a numeric vector, \code{randomEffectsNodes} will be split by
#'         \code{split}(\code{randomEffectsNodes}, \code{control$split}). The
#'         last option allows arbitrary control over how
#'         \code{randomEffectsNodes} are blocked.
#'
#'   \item \code{warn}. If TRUE (default), a warning is issued if
#'   \code{randomEffectsNodes} and/or \code{calcNodes} are provided but have
#'   other or missing elements relative to their defaults.
#'
#'   \item \code{innerOptimControl}. See \code{optimControl}.
#'
#'   \item \code{innerOptimMethod}. See \code{optimMethod}.
#'
#'   \item \code{innerOptimStart}. see \code{optimStart}.
#'
#'   \item \code{outOptimControl}. A list of control parameters for maximizing
#'         the Laplace log-likelihood using \code{optim}. See 'Details' of
#'         \code{\link{optim}} for further information.
#'
#' }
#'
#' @author Wei Zhang, Perry de Valpine
#' 
#' @name laplace
#' 
#' @aliases Laplace buildLaplace
#'
#' @examples 
#' pumpCode <- nimbleCode({ 
#'   for (i in 1:N){
#'     theta[i] ~ dgamma(alpha, beta)
#'     lambda[i] <- theta[i] * t[i]
#'     x[i] ~ dpois(lambda[i])
#'   }
#'   alpha ~ dexp(1.0)
#'   beta ~ dgamma(0.1, 1.0)
#' })
#' pumpConsts <- list(N = 10, t = c(94.3, 15.7, 62.9, 126, 5.24, 31.4, 1.05, 1.05, 2.1, 10.5))
#' pumpData <- list(x = c(5, 1, 5, 14, 3, 19, 1, 1, 4, 22))
#' pumpInits <- list(alpha = 0.1, beta = 0.1, theta = rep(0.1, pumpConsts$N))
#' pump <- nimbleModel(code = pumpCode, name = "pump", constants = pumpConsts, 
#'                     data = pumpData, inits = pumpInits, buildDerivs = TRUE)
#'                     
#' # Build Laplace approximation
#' pumpLaplace <- buildLaplace(pump)
#' 
#' \dontrun{
#' # Compile the model
#' Cpump <- compileNimble(pump)
#' CpumpLaplace <- compileNimble(pumpLaplace, project = pump)
#' # Calculate MLEs of parameters
#' MLEres <- CpumpLaplace$findMLE()
#' # Calculate estimates and standard errors for parameters and random effects on original scale
#' allres <- CpumpLaplace$summary(MLEres, calcRandomEffectsStdError = TRUE)
#' }
#'
#' @references
#'
#' Kass, R. and Steffey, D. (1989). Approximate Bayesian inference in
#' conditionally independent hierarchical models (parametric empirical Bayes
#' models). \emph{Journal of the American Statistical Association}, 84(407),
#' 717–726.
#' 
#' Skaug, H. and Fournier, D. (2006). Automatic approximation of the marginal
#' likelihood in non-Gaussian hierarchical models. \emph{Computational
#' Statistics & Data Analysis}, 56, 699–709.
#' 
NULL
