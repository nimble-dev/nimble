## NIMBLE Laplace approximation
## Laplace base class
APPROX_BASE <- nimbleFunctionVirtual(
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
    },
		get_inner_mode = function(){ returnType(double(1))},
		get_inner_negHessian = function(){returnType(double(2))},
		get_inner_negHessian_chol = function(){returnType(double(2))}
  )
)

## A single Laplace approximation for only one scalar random effect node
buildOneLaplace1D_inner <- function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
  buildOneAGHQuad1D_inner(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart)
}

buildOneAGHQuad1D_inner <- nimbleFunction(
  contains = APPROX_BASE,
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

		## Cache values for posterior approximations via simulation.
		## Mode is called max_inner_logLik_saved_par
		saved_inner_negHess <- matrix(0, nrow = 1, ncol = 1)
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
			saved_inner_negHess <<- matrix(negHessValue, ncol = 1, nrow = 1)
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
			saved_inner_negHess <<- matrix(exp(logdetNegHessian), nrow = 1, ncol = 1)

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
			saved_inner_negHess <<- matrix(exp(logdetNegHessian), nrow = 1, ncol = 1)

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
			saved_inner_negHess <<- matrix(negHessian, nrow = 1, ncol = 1)

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
			saved_inner_negHess <<- matrix(negHessian, nrow = 1, ncol =  1)

      ## invNegHessian <- inverse(negHessian)
      grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p_internal(p, reTransform)
      grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re_internal(p, reTransform)[1]
      hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform)[,1]
      ans <- gr_joint_logLik_wrt_p_internal(p, reTransform) - 
        0.5 * (grlogdetNegHesswrtp + hesslogLikwrtpre * (grlogdetNegHesswrtre / negHessian))
      return(ans)
      returnType(double(1))
    },
		get_inner_mode = function(){ returnType(double(1)); return(max_inner_logLik_saved_par)},
		get_inner_negHessian = function(){ returnType(double(2)); return(saved_inner_negHess)},
		get_inner_negHessian_chol = function(){ returnType(double(2)); return(sqrt(saved_inner_negHess))}
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
buildOneLaplace_inner <- function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
  buildOneAGHQuad_inner(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart)
}

buildOneAGHQuad_inner <- nimbleFunction(
  contains = APPROX_BASE,
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

		## Cache values for access in INLA like methods.
		saved_inner_negHess <- matrix(0, nrow = nre, ncol = nre)
		saved_inner_negHess_chol <- matrix(0, nrow = nre, ncol = nre)
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
    ## Logdet negative Hessian
    logdetNegHess_internal = function(p = double(1), reTransform = double(1)) {
      negHessian <- negHess(p, reTransform)
			saved_inner_negHess <<- negHessian
      cholNegHess <- chol(negHessian)
			saved_inner_negHess_chol <<- cholNegHess
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
			saved_inner_negHess_chol <<- chol_negHess
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
      logdetNegHessian <- logdetNegHess_internal(p, reTransform)
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
      logdetNegHessian <- logdetNegHess_internal(p, reTransform)
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
    },
		get_inner_mode = function(){ returnType(double(1)); return(max_inner_logLik_saved_par)},
		get_inner_negHessian = function(){ returnType(double(2)); return(saved_inner_negHess)},
		get_inner_negHessian_chol = function(){ returnType(double(2)); return(saved_inner_negHess_chol)}
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


## Inner quadrature for posterior inference.
buildApproxPosterior <- nimbleFunction(
  name = 'AGHQuad',
  setup = function(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes, calcNodesOther,
                   control = list(), approxMethods = list(hyperGrid = 'ccd', nQuadAGHQ = 5)) {
    if(is.null(control$split)) split <- TRUE else split <- control$split
    if(is.null(control$check))   check <- TRUE else  check <- control$check
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
                                  check = check)
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
    optimControlArgNames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
	
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
    AGHQuad_nfl <- nimbleFunctionList(APPROX_BASE)
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
        if(nre > 1) {
					if(nQuad > 1) nimCat('  [Note] Adaptive Gauss-Hermite quadrature is not implemented for      multivariate integration.\n Defaulting to Laplace approximation.')
					AGHQuad_nfl[[1]] <- buildOneAGHQuad_inner(model, paramNodes, randomEffectsNodes, calcNodes,
							innerOptControl, innerOptMethod, innerOptStart)
        } else { 
					AGHQuad_nfl[[1]] <-  buildOneAGHQuad1D_inner(model, paramNodes, randomEffectsNodes, calcNodes, innerOptControl, "CG", innerOptStart)
				}  
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
						if(nQuad > 1) nimCat('  [Note] Adaptive Gauss-Hermite quadrature is not implemented for multivariate integration.\n Defaulting to Laplace Approximation.')

								AGHQuad_nfl[[i]] <-  buildOneAGHQuad_inner(model, paramNodes, these_reNodes, these_calcNodes, innerOptControl, innerOptMethod, innerOptStart)
							}
							else AGHQuad_nfl[[i]] <-  buildOneAGHQuad1D_inner(model, paramNodes, these_reNodes, these_calcNodes, innerOptControl, "CG", innerOptStart)
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
    
    paramNodesAsScalars <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
    npar <- length(paramNodesAsScalars)
    paramNodesAsScalars_vec <- paramNodesAsScalars
    if(npar == 1) paramNodesAsScalars_vec <- c(paramNodesAsScalars, "_EXTRA_")
    paramNodesAsScalars_first <- paramNodesAsScalars[1]
    if(npar == 1) p_indices <- c(1, -1)
    else p_indices <- 1:npar
    ## setupOutputs(reNodesAsScalars, paramNodesAsScalars)
    
    ## Automated transformation for parameters
    paramsTransform <- parameterTransformI(model, paramNodes, control = list(allowDeterm = allowNonPriors))
    pTransform_length <- paramsTransform$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
   
		## Define the hyperparameter grid:
		## Probably need to say if(pTransform_length > 0) but surely that is implied...
		I_CCD <- 1
		I_AGHQ <- 2
		## To be able to compare accuracy we will build both grids, or at least an empty version of aghq if not requested.
    theta_grid_nfl <- nimbleFunctionList(GRID_BASE)
		if(approxMethods$hyperGrid == 'ccd'){
			gridMethod <- I_CCD
			theta_grid_nfl[[I_CCD]] <- buildQuadGrid(d = pTransform_length, nQ = 1, method = I_CCD, nre = nre)
			theta_grid_nfl[[I_AGHQ]] <- buildQuadGrid(d = pTransform_length, nQ = 1, method = I_AGHQ, nre = nre)
		}else{
			gridMethod <- I_AGHQ
			theta_grid_nfl[[I_AGHQ]] <- buildQuadGrid(d = pTransform_length, nQ = approxMethods$nQuadAGHQ, method = I_AGHQ, nre = nre)					
		}
		gridBuilt <- c(0,0)
		thetaModeIndex <- 1	## Default for CCD. Will update when building grid.
		nQuadAGHQ <- approxMethods$nQuadAGHQ
		nGrid <- theta_grid_nfl[[gridMethod]]$getGridSize()
		wgtA <- 0

		# zGrid <- gridPts$z_nodes
		# zWeights <- gridPts$weights
		# nGrid <- nrow(zGrid)
		# theta_grid_nfl[[gridMethod]] <- storeGridValues(z = zGrid, wgt = zWeights, nre = nre)

    ## Cache values from optim step.
		pTransformPostMode <- rep(0, pTransform_length)
    logPostProbMode <- 0
		negHesspTransformPostMode <- matrix(0,  nrow = pTransform_length, ncol = pTransform_length)
    
		## Cache values for marginal distributions in AGHQ over d-1.
		pTransformFix <- 0
		indexFix <- 0
		
		## Need to have eigen decomp saved.
		hessEigenVecs <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		hessEigenVals <- numeric(pTransform_length)
		ATransform <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		AInverse <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		covTheta <- matrix(0, nrow=pTransform_length, ncol=pTransform_length)
		
		## We will want to cache the standard deviation skew terms.
		skewedStdDev <- matrix(0, nrow = pTransform_length, ncol = 2)
		extraPoints <- c(-15.0, -10.0, -7.0, 7.0, 10.0, 15.0, -5.0, -3.0, -2.0, -1.0, -0.5, -0.25, 0.0, 0.25, 0.5, 1.0, 2.0, 3.0, 5.0)
		aghqPoints <- pracma::gaussHermite(51)$x*sqrt(2)
		zMargGrid <- sort(unique(c(extra.points, aghq.points)))

		## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    ## Default calculation method for AGHQuad
    methodID <- 2
    ## The nimbleList definitions AGHQuad_params and AGHQuad_summary
    ## have moved to predefined nimbleLists.
  },## End of setup
  run = function(){},
  methods = list(
    getNodeNamesVec = function(returnParams = logical(0, default = TRUE)) {
      returnType(character(1))
      if(returnParams) return(paramNodesAsScalars_vec)
      else return(reNodesAsScalars_vec)
    },
    getNodeNameSingle = function(returnParams = logical(0, default = TRUE)) {
      returnType(character())
      if(returnParams) return(paramNodesAsScalars_first)
      else return(reNodesAsScalars_first)
    },
    setMethod = function(method = integer()) {
      if(nre == 0) print("AGHQuad or Laplace approximation is not needed for the given model: no random effects")
      if(!any(c(1, 2, 3) == method)) stop("Choose a valid method ID from 1, 2, and 3")
      methodID <<- method
    },
    getMethod = function() {
      return(methodID)
      returnType(integer())
    },
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(pTransform_length == 1){
        if(length(pTransform_indices) == 2){
          pTransform_indices <<- numeric(length = 1, value = 1)
		  pTransformPostMode <<- numeric(length = 1, value = 0)
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
		## Added functions for INLA Grid
		##***************************************	
		buildHyperGrid = function(){
			theta_grid_nfl[[gridMethod]]$buildGrid()
			gridBuilt[gridMethod] <<- 1
			thetaModeIndex <<- theta_grid_nfl[[gridMethod]]$getThetaModeIndex()
			wgtA <<- 0
		},
		changeHyperGrid = function(method = character(0, default = "aghq"), 
				nQUpdate = integer(0, default = 3)){
			if(method == "ccd"){
				gridMethod <<- I_CCD
				if(gridBuilt[gridMethod] == 0) {
					theta_grid_nfl[[gridMethod]]$buildGrid()
					wgtA <<- 0
				}	
			} else {
				gridMethod <<- I_AGHQ
				if(gridBuilt[gridMethod] == 1 & nQUpdate == nQuadAGHQ){
					## Do nothing.
				}else{
					theta_grid_nfl[[gridMethod]]$resetGrid(nQUpdate = nQUpdate)	## Only use resetGrid for aghq.
					nQuadAGHQ <<- nQUpdate
					wgtA <<- 0
				}
			}
			gridBuilt[gridMethod] <<- 1
			nGrid <<- theta_grid_nfl[[gridMethod]]$getGridSize()
			thetaModeIndex <<- theta_grid_nfl[[gridMethod]]$getThetaModeIndex()			
		},
		## Added functions for INLA to get posterior density
		##***************************************			
		calcPrior_p = function(p = double(1)){
			values(model, paramNodes) <<- p
			ans <- model$calculate(paramNodes)	## Add updates to deterministic nodes!
			return(ans)
			returnType(double())
		},
		calcPrior_pTransformed = function(pTransform = double(1)) {
			p <- paramsTransform$inverseTransform(pTransform)
			values(model, paramNodes) <<- p
			ans <- model$calculate(paramNodes)
			return(ans)
			returnType(double())
		},
		calcPostLogProb = function(p = double(1), trans = logical(0, default = FALSE)) {
			if(trans){
				pstar <- paramsTransform$inverseTransform(p)
			}else{ 
				pstar <- p
			}
			ans <- calcLogLik(pstar) + calcPrior_p(pstar)
			return(ans)
			returnType(double())
		},
		calcPostLogProb_pTransformed = function(pTransform = double(1)) {
			ans <- calcPostLogProb(pTransform, trans = TRUE) + logDetJacobian(pTransform)
			returnType(double())
			return(ans)
		},
		## Use internally to loop through quadrature point evaluations.
		calcPostLogProb_pTransformed_multi = function(pTransform = double(2)) {
			dm <- dim(pTransform)[1]
			ans <- numeric(length = dm)
			for( i in 1:dm){
				ans[i] <- calcPostLogProb_pTransformed(pTransform[i,])
			}
			returnType(double(1))
			return(ans)
		},
		calcPostLogProb_pTransformedj = function(pTransform = double(1))
		{
			vals <- pTransformPostMode
			for( i in 1:pTransform_length){
				if( i < indexFix ) vals[i] <- pTransform[i]			
				if( i == indexFix ) vals[i] <- pTransformFix
				if( i > indexFix ) vals[i] <- pTransform[i-1]
			}
			ans <- calcPostLogProb_pTransformed(vals)	
			returnType(double())
			return(ans)
		},
		logDetJacobian = function(pTransform = double(1)){
			ans <- paramsTransform$logDetJacobian(pTransform)
			return(ans)
			returnType(double())
		},
		##-----------------------------------------------------
		## Derivatives:
		## Transform parameters to real scale
		gr_logDetJacobian = function(pTransform = double(1))
		{
			ans <- derivs(logDetJacobian(pTransform), wrt = pTransform_indices, order = 1)
			return(ans$jacobian[1,])
			returnType(double(1))
		},
		gr_prior = function(p = double(1))
		{
			ans <- derivs(calcPrior_p(p), wrt = p_indices, order = 1)
			return(ans$jacobian[1,])
			returnType(double(1))
		},
		gr_postLogProb_pTransformed = function(pTransform = double(1))
		{
			## *** Repeated gradients and inverse.
			pDerivs <- derivs_pInverseTransform(pTransform, c(0, 1))
			grLogDetJacobian <- gr_logDetJacobian(pTransform)
			grLogLikTrans <- gr_logLik(pTransform, TRUE)

			p <- pDerivs$value
			grPrior <- gr_prior(p)
			grPriorTrans <- (grPrior %*% pDerivs$jacobian)[1,]
			
			ans <- grLogLikTrans + grPriorTrans + grLogDetJacobian
			return(ans)
			returnType(double(1))
		},
		gr_postLogProb_pTransformedj = function(pTransform = double(1))
		{
			vals <- pTransformPostMode
			for( i in 1:pTransform_length){
				if( i < indexFix ) vals[i] <- pTransform[i]			
				if( i == indexFix ) vals[i] <- pTransformFix
				if( i > indexFix ) vals[i] <- pTransform[i-1]
			}
			ans <- gr_postLogProb_pTransformed(vals)
			ansj <- ans[pTransform_indices != indexFix]
			return(ansj)
			returnType(double(1))
		},
		spectralTransformation = function(negHess = double(2)){
			## Save Eigen Vectors and Values.
			eigenDecomp <- nimEigen(negHess)
			hessEigenVecs <<- eigenDecomp$vectors
			hessEigenVals <<- eigenDecomp$values		
			ATransform <<- hessEigenVecs %*% diag(1/sqrt(hessEigenVals))
			AInverse <<- diag(sqrt(hessEigenVals)) %*% t(hessEigenVecs)
			covTheta <<- ATransform %*% ATransform	## Efficient to do this here.
		},
		## Transform z to theta.
		z_to_theta = function(z = double(1)) {
			theta <- pTransformPostMode + (ATransform %*% z)
			returnType(double(1))
			return(theta[,1])
		},		
		## Functions to approximate posterior distribution of theta.
		theta_to_z = function(theta = double(1)) {
			z <- AInverse %*% (theta - pTransformPostMode)
			returnType(double(1))
			return(z[,1])
		},				
		## Functions to approximate posterior distribution of theta.
		calcSkewedSD = function() {
			for( i in 1:pTransform_length){
				z <- numeric(value = 0, length = pTransform_length)
				z[i] <- -sqrt(2)
				theta <- z_to_theta(z)
				logDens2Neg <- calcPostLogProb_pTransformed(theta)
				skewedStdDev[i, 1] <<- sqrt(2 / (2.0 * (logPostProbMode-logDens2Neg))) 	## numerator (-sqrt(2)) ^2
				z[i] <- sqrt(2)
				theta <- z_to_theta(z)
				logDens2Pos <- calcPostLogProb_pTransformed(theta)
				skewedStdDev[i, 2] <<- sqrt(2 / (2.0 * (logPostProbMode-logDens2Pos))) 	## numerator (-sqrt(2)) ^2
			}
		},
		calcMarginalLogLikApprox = function(){
			logDetNegHess <- sum(log(hessEigenVals)) ## 1/det|H|^0.5.		
			marg <- logPostProbMode + 0.5*pTransform_length*log(2*pi) - 
				0.5 * logDetNegHess + sum(log((skewedStdDev[,1] + skewedStdDev[,2]))/2)
			returnType(double())
			return(marg)
			## Line 2748 in r-inla/blob/devel/gmrflib/approx-inference.c
			## Commit # ef4eb20
			# + sum( log((SD[,1] * SD[,2])) )	## Version from INLA but doesn't make sense to me.
		},
		## This is a loop for calculating theta on the grid points.
		## It stores all values we need for simulation inference on the fixed and random-effects.
		## Once this is called, we can then simulate.
		calcHyperGrid = function(){
			if( gridBuilt[gridMethod] == 0 ) buildHyperGrid()
			for( i in 1:nGrid ){
				if(theta_grid_nfl[[gridMethod]]$calcCheck(i) == 0 ){
					z <- theta_grid_nfl[[gridMethod]]$getNodes(i)
					theta <- z_to_theta(z)
					# theta <- theta_grid_nfl[[gridMethod]]$getTheta(i)
					logDensi <- calcPostLogProb_pTransformed(theta)
					tmpChol <- get_inner_cholesky()								# Get the updated vals from Laplace.
					tmpMode <- get_inner_mode()
					theta_grid_nfl[[gridMethod]]$saveValues(i=i, theta = theta, 
											logDensity = logDensi,
											innerCholesky = tmpChol, innerMode = tmpMode)
				}
			}
		},
		## Quadrature based marginal log-likelihood
		## Unlikely to be particularly accurate for CCD.
		## Make sure to add 1/det|H|^0.5.
		calcMarginalLogLikQuad = function(){
			marg <- 0
			logDetNegHess <- sum(log(hessEigenVals)) ## Adjust weights by 1/det|H|^0.5.
			for(i in 1:nGrid){
				lik <- exp(theta_grid_nfl[[gridMethod]]$getLogDensity(i) - logPostProbMode)
				marg <- marg + lik*theta_grid_nfl[[gridMethod]]$getWeights(i)
			}
			returnType(double())
			return(log(marg) + logPostProbMode - 0.5*logDetNegHess)
		},
		## Use this to test and inspect. Might keep.
		getQuadratureGrid = function(){
			nGrid <<- theta_grid_nfl[[gridMethod]]$getGridSize()
			ans <- matrix(0, nrow = nGrid, ncol = pTransform_length + 1 + 1)
			for(i in 1:nGrid){
				ans[i, 1:pTransform_length] <- theta_grid_nfl[[gridMethod]]$getTheta(i)
				ans[i, pTransform_length + 1] <- theta_grid_nfl[[gridMethod]]$getWeights(i)
				ans[i, pTransform_length + 2] <- theta_grid_nfl[[gridMethod]]$getLogDensity(i)
			}
			returnType(double(2))
			return(ans)
		},
		findMarginalHyperAGHQ = function(pIndex = integer(), size = integer())
		{
		## Some efficiencies here if nQuad is the same as the previously built grid.
		## if(nAGHQ == theta_grid_nfl[[I_AGHQ]$getGridSize())
		## Need to build a d - 1 AGHQ grid object in setup.
		## findPostModeFixedj()
		## Then calculate the log density
		},
		findMarginalHyperIntFree = function(pIndex = integer())
		{
			std_dev <- sqrt(covTheta[pIndex, pIndex])
			thetai <- rep(0, pTransform_length)
			logDens <- rep(0, 69)
			for( i in 1:69 ){	# Know fixed # of points
				thetai[pIndex] <- pTransformPostMode[pIndex] + zMargGrid[i] * std_dev
				## Find the conditional mean:
				for( j in 1:pTransform_length ){
					if(j != pIndex){
						thetai[j] <- pTransformPostMode[j] + covTheta[pIndex, j]/covTheta[pIndex, pIndex] * 
							(thetai[pIndex] - pTransformPostMode[pIndex])
					}
				}
				## Calculate asymmetric Gaussian:
				zi <- theta_to_z(thetai)
				for( j in 1:pTransform_length ){
					side <- 2
					if(zi[j] <= 0) side <- 1
					logDens[i] <- logDens[i] - 0.5*zi[j]^2/skewedStdDev[j, side]
				}
			}
			## Note to self: Explore normalization... Should we use the spline? Can we know it analytically?
			returnType(double(1))
			return(logDens)			
		},		
		findApproxPosterior = function(){
		
			## Basic approx posterior STEPS:
			##-------------------------------------
			## 1) Build the hyperparameter grid.
			if( gridBuilt[gridMethod] == 0 ) buildHyperGrid()
			nGrid <<- theta_grid_nfl[[gridMethod]]$getGridSize()

			## 2) Find posterior mode: 
			## Values are saved to theta_grid_nfl and locally cached.
			## Should think about initialization.
			findPostMode(pStartTransform = rep(0, pTransform_length), method = "BFGS", hessian = TRUE)

			## 3) Calculate skew and if gridMethod == CCD, skew grid.
			calcSkewedSD()
			if(gridMethod == I_CCD) theta_grid_nfl[[I_CCD]]$skewGridPoints(skewedStdDev)

			## 4) Calculate the density on the grid points. Saves to theta_grid_nfl
			## This is used for inference on the fixed and random-effects.
			calcHyperGrid()
			
			## 5) Calculate Marginal log-Likelihood
			## Based on asymetric Gaussian assumption of the marginal of theta.
			marginalAG <- calcMarginalLogLikApprox()
			marginalQuad <- calcMarginalLogLikQuad()	## This one probably make sense only for AGHQ

			## 6) Marginals for theta: 
			## Options are integration free or AGHQ Stringer Method.
			## Need to do marginals aghq for stringer and the integration free method.
			## Nimble R Call:
			
			## 7) Marginals for Fixed and Random-Effects: 
			## Only simulation based.
			## Need to simulate based on the quadrature grid.

		},
		##***************************************
		## Outer Optim (MAP) for INLA:
		##**************************************
		findPostMode = function(pStartTransform  = double(1, default = Inf),
											 method  = character(0, default = "BFGS"),
											 hessian = logical(0, default = TRUE)) {
			if(any(abs(pStartTransform) == Inf)){
				pStart <- values(model, paramNodes)
				pStartTransform <- paramsTransform$transform(pStart)
			}
			## In case bad start values are provided 
			if(any_na(pStartTransform) | any_nan(pStartTransform) | any(abs(pStartTransform)==Inf)) pStartTransform <- rep(0, pTransform_length)
			optRes <- optim(pStartTransform, calcPostLogProb_pTransformed, gr_postLogProb_pTransformed, method = method, 
				control = outOptControl, hessian = hessian)
			if(optRes$convergence != 0) 
				print("Warning: optim has a non-zero convergence code: ", optRes$convergence, ".\n",
							"The control parameters of optim can be adjusted in the control argument of\n",
							"buildLaplace or buildAGHQuad via list(outOptimControl = list()).")
			## Back transform results to original scale
			# optRes$par <- paramsTransform$inverseTransform(optRes$par)
			pTransformPostMode <<- optRes$par
			logPostProbMode <<- optRes$value
			negHesspTransformPostMode <<- -optRes$hessian

			## Store these points in the quadrature object for use later.
			## This does the mode only.
			tmpChol <- get_inner_cholesky()
			tmpMode <- get_inner_mode()

			if( gridBuilt[gridMethod] == 0 ) buildHyperGrid()
			## i = 0 defaults to mode.
			theta_grid_nfl[[gridMethod]]$saveValues(i=0, theta = pTransformPostMode, 
									logDensity = logPostProbMode,
									innerCholesky = tmpChol, innerMode = tmpMode)

			## Save all the spectral transformation information.
			spectralTransformation(negHesspTransformPostMode)

			return(optRes)
			returnType(optimResultNimbleList())
		},
		findPostModeFixedj = function(pStartTransform  = double(1, default = Inf),
						 j = integer(0, default = 1), 
						 pTransformFixed = double(0, default = 0),
											 method  = character(0, default = "BFGS"),
											 hessian = logical(0, default = TRUE)) {
			indexFix <<- j
			pTransformFix <<- pTransformFixed
		
			if(any(abs(pStartTransform) == Inf)){
				pStartTransform <- pTransformPostMode[pTransform_indices != j]
			}
			optRes <- optim(pStartTransform, calcPostLogProb_pTransformedj, gr_postLogProb_pTransformedj, method = method,
					control = outOptControl, hessian = hessian)
			if(optRes$convergence != 0) 
				print("Warning: optim has a non-zero convergence code: ", optRes$convergence, ".\n",
							"The control parameters of optim can be adjusted in the control argument of\n",
							"buildLaplace or buildAGHQuad via list(outOptimControl = list()).")
			
			return(optRes)
			returnType(optimResultNimbleList())
		},	
		##-----------------------------------------------------	
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
		## For external access:
		p_Transform = function(p = double(1))
		{
			pTransform <- paramsTransform$transform(p)
			return(pTransform)
			returnType(double(1))
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
			# stop("Wrong length for pStart in findMLE.")
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
    ## Grab the inner Cholesky from the cached last values.
    get_inner_cholesky = function(){
      if(nre == 0) stop("No random effects in the model")
      cholesky <- matrix(value = 0, nrow = nre, ncol = nre)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        numre <- lenInternalRENodeSets[i]
        cholesky[(tot+1):(tot+numre), (tot+1):(tot+numre)] <- AGHQuad_nfl[[i]]$get_inner_negHessian_chol()
        tot <- tot + numre
      }
      return(cholesky)
      returnType(double(2))
    },
    get_inner_mode = function(){
      if(nre == 0) stop("No random effects in the model")
      raneff <- numeric(nre)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        numre <- lenInternalRENodeSets[i]
        raneff[(tot+1):(tot+numre)] <- AGHQuad_nfl[[i]]$get_inner_mode()
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
    summary = function(MLEoutput                 = optimResultNimbleList(),
                       originalScale             = logical(0, default = TRUE),
                       randomEffectsStdError = logical(0, default = FALSE),
                       jointCovariance       = logical(0, default = FALSE)){
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
          if(jointCovariance) {
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
          if(jointCovariance) ans$vcov <- vcov_pTransform
          else ans$vcov <- matrix(0, nrow = 0, ncol = 0)
        }
      }
      else{
        ## Random effects
        optreTransform <- optimRandomEffects(pTransform)
        optre <- reInverseTransform(optreTransform)
        ntot <- npar + nre
        if(jointCovariance) {
          ## Inverse of the negative Hessian of log-likelihood wrt transformed random effects at MLEs
          inv_negHess <- inverse_negHess(p, optreTransform)
          jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
          #jointInvNegHessZero[1:nre, 1:nre] <- inv_negHess
          jointInvNegHessZero[(npar+1):ntot, (npar+1):ntot] <- inv_negHess
          ## Hessian of log-likelihood wrt to params and transformed random effects
          hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
          ## Derivative of inverse transformation for params
          derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
          JacobpInvTransform   <- derivspInvTransform$jacobian
          ## Jacobian of optimized random effects wrt transformed parameters
          JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
          jointJacob <- matrix(init = FALSE, nrow = ntot, ncol = npar)
          #jointJacob[1:nre, 1:npar] <- JacobOptreWrtParams
          jointJacob[(npar+1):ntot, 1:npar] <- JacobOptreWrtParams
          #jointJacob[(nre+1):ntot, 1:npar] <- diag(npar)
          jointJacob[1:npar, 1:npar] <- diag(npar)
          ## Joint covariance matrix on transformed scale
          vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
          if(originalScale){
            derivs_reInvTransform <- derivs_reInverseTransform(optreTransform, c(0, 1))
            Jacob_reInvTransform  <- derivs_reInvTransform$jacobian
            Jacob_JointInvTransform <- matrix(0, nrow = ntot, ncol = ntot)
            #Jacob_JointInvTransform[1:nre, 1:nre] <- Jacob_reInvTransform
            Jacob_JointInvTransform[(npar+1):ntot, (npar+1):ntot] <- Jacob_reInvTransform
            #Jacob_JointInvTransform[(nre+1):ntot, (nre+1):ntot] <- JacobpInvTransform
            Jacob_JointInvTransform[1:npar, 1:npar] <- JacobpInvTransform
            vcov <- Jacob_JointInvTransform %*% vcov_Transform %*% t(Jacob_JointInvTransform)
            stdErr_p_re <- sqrt(diag(vcov))
            stdErr_p <- stdErr_p_re[1:npar]
            if(randomEffectsStdError){
              ranres$stdErrors <- stdErr_p_re[(npar+1):ntot]
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
            if(randomEffectsStdError){
              stdErr_reTransform <- sqrt(diag(vcov_Transform)[(npar+1):ntot])
              ranres$stdErrors <- stdErr_reTransform
            }
            else{
              ranres$stdErrors <- numeric(0)
            }
            ans$vcov <- vcov_Transform
            pres$estimates <- pTransform
            pres$stdErrors <- sqrt(diag(vcov_Transform)[1:npar])
            ranres$estimates <- optreTransform
          }
        }## End of if(jointCovariance)
        else { ## Do not return joint covariance matrix
          if(originalScale){## On original scale
            pres$estimates <- p
            ranres$estimates <- optre
            if(randomEffectsStdError){
              ## Joint covariance matrix on transform scale
              inv_negHess <- inverse_negHess(p, optreTransform)
              # jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
              # jointInvNegHessZero[1:nre, 1:nre] <- inv_negHess
              ## Hessian of log-likelihood wrt to params and transformed random effects
              hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
              ## Derivative of inverse transformation for params
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              ## Covariance matrix for params on the original scale
              vcov_p <- JacobpInvTransform %*% vcov_pTransform %*% t(JacobpInvTransform)
              ## Jacobian of optimized random effects wrt transformed parameters
              JacobOptreWrtParams <- inv_negHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
              # jointJacob <- matrix(NA, nrow = ntot, ncol = npar)
              # jointJacob[1:nre, 1:npar] <- JacobOptreWrtParams
              # jointJacob[(nre+1):ntot, 1:npar] <- diag(npar)
              ## Join covariance matrix on transformed scale
              # vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
              ## Covariance matrix for random effects (transformed) 
              vcov_reTransform <- inv_negHess + JacobOptreWrtParams %*% vcov_pTransform %*% t(JacobOptreWrtParams)
              ## Derivatives information
              derivs_reInvTransform <- derivs_reInverseTransform(optreTransform, c(0, 1))
              Jacob_reInvTransform  <- derivs_reInvTransform$jacobian
              # Jacob_JointInvTransform <- matrix(0, nrow = ntot, ncol = ntot)
              # Jacob_JointInvTransform[1:nre, 1:nre] <- Jacob_reInvTransform
              # Jacob_JointInvTransform[(nre+1):ntot, (nre+1):ntot] <- JacobpInvTransform
              stdErr_re <- numeric(nre)
              for(i in 1:nre){
                var_i <- (Jacob_reInvTransform[i,,drop=FALSE] %*% vcov_reTransform %*% t(Jacob_reInvTransform[i,,drop=FALSE]))[1,1]
                stdErr_re[i] <- sqrt(var_i)
              }
              stdErr_p <- sqrt(diag(vcov_p))
              pres$stdErrors   <- stdErr_p
              ranres$stdErrors <- stdErr_re
              ans$vcov <- vcov_p
            }## End of if(randomEffectsStdError)
            else { ## Do not calculate standard errors of random effects estimates
              derivspInvTransform  <- derivs_pInverseTransform(pTransform, c(0, 1))
              JacobpInvTransform   <- derivspInvTransform$jacobian
              ## Covariance matrix for params on the original scale
              vcov_p <- JacobpInvTransform %*% vcov_pTransform %*% t(JacobpInvTransform)
              # stdErr_p <- numeric(npar)
              # for(i in 1:npar){
              #   var_p_i <- (JacobpInvTransform[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobpInvTransform[i,,drop=FALSE]))[1,1]
              #   stdErr_p[i] <- sqrt(var_p_i)
              # }
              stdErr_p <- sqrt(diag(vcov_p))
              pres$stdErrors <- stdErr_p
              ranres$stdErrors <- numeric(0)
              ans$vcov <- vcov_p
            }
          }## End of if(originalScale)
          else {## On transformed scale
            pres$estimates <- pTransform
            pres$stdErrors <- stdErr_pTransform
            ranres$estimates <- optreTransform
            ans$vcov <- vcov_pTransform
            if(randomEffectsStdError){
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
                     gr_otherLogLik_internal = list(),
										 logDetJacobian = list(),
										 calcPrior_p = list()
										)
)