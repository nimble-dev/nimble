## NIMBLE Laplace approximation
## Laplace base class
#' @rdname laplace
#' @export
Laplace_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
    Laplace1 = function(p = double(1)){
      returnType(double())
    },
    Laplace2 = function(p = double(1)){
      returnType(double())
    },
    Laplace3 = function(p = double(1)){
      returnType(double())
    },
    gr_Laplace1 = function(p = double(1)){
      returnType(double(1))
    },
    gr_Laplace2 = function(p = double(1)){
      returnType(double(1))
    },
    gr_Laplace3 = function(p = double(1)){
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
#' @rdname laplace
#' @export
nimOneLaplace1D <- nimbleFunction(
  contains = Laplace_BASE,
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
    ## Check the number of random effects is 1
    nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))
    if(length(nre) != 1) stop("Number of random effects for nimOneLaplace1D must be 1.")
    ## Check and add necessary deterministic nodes into calcNodes
    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
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
    inner_derivsInfo    <- nimble:::makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- nimble:::makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes
    
    ## Automated transformation for random effects to ensure range of (-Inf, Inf) 
    reTrans <- parameterTransform(model, randomEffectsNodes)
    
    ## The following are used for caching values and gradient in the Laplace3 system
    Laplace3_saved_value <- numeric(1)
    Laplace3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    Laplace3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    ## The following are used for caching values and gradient in the Laplace3 system
    max_inner_logLik_saved_par <- as.numeric(c(1, -1))
    max_inner_logLik_saved_value <- numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    ## Record the maximum Laplace loglikelihood value for obtaining inner optimization start values
    max_Laplace <- -Inf 
    max_Laplace_saved_re_value <- as.numeric(c(1, -1))
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
      else if(startID == 2) ans <- max_Laplace_saved_re_value
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
      max_Laplace_saved_re_value <<- fix_one_vec(max_Laplace_saved_re_value)
      if(startID == 3) optStart <<- fix_one_vec(optStart)
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        Laplace3_saved_gr <<- fix_one_vec(Laplace3_saved_gr)
        Laplace3_previous_p <<- fix_one_vec(Laplace3_previous_p)
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
        print("Warning: optim does not converge for the inner optimization of Laplace approximation")
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
        print("Warning: optim does not converge for the inner optimization of Laplace approximation")
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
    update_Laplace3_with_gr = function(p = double(1), reset = logical(0, default = FALSE)) {
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
      
      Laplace_value <- maxValue - 0.5 * logdetNegHess_value + 0.5 * 1 * log(2*pi)
      Laplace3_saved_value <<- Laplace_value
      
      gr_Laplace_v <- gr_logLik_wrt_p - 0.5*(gr_logdetNegHess_wrt_p_v + hess_cross_terms * (gr_logdetNegHess_wrt_re_v / negHessValue))
      Laplace3_saved_gr <<- gr_Laplace_v
      return(ans$value)
      returnType(double(1))
    },
    Laplace3_update = function(p = double(1)) {
      if(any(p != Laplace3_previous_p)) {
        update_Laplace3_with_gr(p)
        Laplace3_previous_p <<- p
      }
    },
    Laplace3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      Laplace3_update(p)
      if(Laplace3_saved_value > max_Laplace) {
        max_Laplace <<- Laplace3_saved_value
        max_Laplace_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(Laplace3_saved_value)
      returnType(double())
    },
    gr_Laplace3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      Laplace3_update(p)
      return(Laplace3_saved_gr)
      returnType(double(1))
    },
    ## Laplace approximation 2: double tapping with separate components
    Laplace2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ## Laplace approximation
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * 1 * log(2*pi)
      if(ans > max_Laplace) {
        max_Laplace <<- ans
        max_Laplace_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Laplace approximation 1: single tapping with separate components
    Laplace1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ## Laplace approximation
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * 1 * log(2*pi)
      if(ans > max_Laplace) {
        max_Laplace <<- ans
        max_Laplace_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation (version 2) w.r.t. parameters
    gr_Laplace2 = function(p = double(1)){
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
    gr_Laplace1 = function(p = double(1)){
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
) ## End of nimOneLaplace1D


## A single Laplace approximation for models with more than one scalar random effect node
#' @rdname laplace
#' @export
nimOneLaplace <- nimbleFunction(
  contains = Laplace_BASE,
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes, optimControl, optimMethod, optimStart) {
    ## Check and add necessary deterministic nodes into calcNodes
    paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
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
    inner_derivsInfo    <- nimble:::makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- nimble:::makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes
    
    ## The following are used for caching values and gradient in the Laplace3 system
    Laplace3_saved_value <- -Inf #numeric(1)
    Laplace3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    Laplace3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    
    max_inner_logLik_saved_par <- if(nreTrans > 1) numeric(nreTrans) else as.numeric(c(1, -1))
    max_inner_logLik_saved_value <- -Inf #numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    
    ## Record the maximum Laplace loglikelihood value for obtaining inner optimization start values
    max_Laplace <- -Inf 
    max_Laplace_saved_re_value <- if(nreTrans > 1) numeric(nreTrans) else as.numeric(c(1, -1))
    
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
      else if(startID == 2) ans <- max_Laplace_saved_re_value
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
        max_Laplace_saved_re_value <<- fix_one_vec(max_Laplace_saved_re_value)
      }
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        Laplace3_saved_gr <<- fix_one_vec(Laplace3_saved_gr)
        Laplace3_previous_p <<- fix_one_vec(Laplace3_previous_p)
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
        print("Warning: optim does not converge for the inner optimization of Laplace approximation")
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
        print("Warning: optim does not converge for the inner optimization of Laplace approximation")
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
    update_Laplace3_with_gr = function(p = double(1), reset = logical(0, default = FALSE)) {
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
      
      Laplace_value <- maxValue - 0.5 * logdetNegHess_value + 0.5 * nreTrans * log(2*pi)
      Laplace3_saved_value <<- Laplace_value
      
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
      gr_Laplace_v <- gr_logLik_wrt_p - 0.5*(gr_logdetNegHess_wrt_p_v + v %*% w )
      # print( gr_Laplace_v )
      Laplace3_saved_gr <<- numeric(gr_Laplace_v, length = npar)
      return(ans$value)
      returnType(double(1))
    },
    Laplace3_update = function(p = double(1)) {
      if(any(p != Laplace3_previous_p)) {
        update_Laplace3_with_gr(p)
        Laplace3_previous_p <<- p
      }
    },
    Laplace3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      Laplace3_update(p)
      ans <- Laplace3_saved_value
      if(ans > max_Laplace) {
        max_Laplace <<- ans
        max_Laplace_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    gr_Laplace3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      Laplace3_update(p)
      return(Laplace3_saved_gr)
      returnType(double(1))
    },
    ## Laplace approximation 2: double tapping with separate components
    Laplace2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      if(ans > max_Laplace) {
        max_Laplace <<- ans
        max_Laplace_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Laplace approximation 1: single tapping with separate components
    Laplace1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed
      logdetNegHessian <- logdetNegHess(p, reTransform)
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      if(ans > max_Laplace) {
        max_Laplace <<- ans
        max_Laplace_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation 2 w.r.t. parameters
    gr_Laplace2 = function(p = double(1)){
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
    gr_Laplace1 = function(p = double(1)){
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
) ## End of nimOneLaplace

## Main function for Laplace approximation
#' @rdname laplace 
#' @export
buildLaplace <- nimbleFunction(
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes, control = list()) {
    if(is.null(control$split)) split <- TRUE else split <- control$split
    if(is.null(control$warn)) warn <- TRUE else warn <- control$warn
    
    paramProvided <- !missing(paramNodes)
    reProvided    <- !missing(randomEffectsNodes)
    calcProvided  <- !missing(calcNodes)
    
    if(!paramProvided) {
      paramNodes <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE)
    } else {
      paramNodes <- model$expandNodeNames(paramNodes)
    }
    
    if((!reProvided) || warn) {
      paramDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
      reNodesDefault <- model$getNodeNames(latentOnly = TRUE)
      reNodesDefault <- intersect(reNodesDefault, paramDeps)
    }
    if(reProvided){
      randomEffectsNodes <- model$expandNodeNames(randomEffectsNodes)
    }
    if(reProvided && warn) {
      reCheck <- setdiff(reNodesDefault, randomEffectsNodes)
      if(length(reCheck)) {
        errorNodes <- paste0(head(reCheck, n = 4), sep = ", ", collapse = ", ")
        if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some random effects (latent states) in the model that look\n",
                       "like they should be included in randomEffectsNodes for Laplace approximation\n",
                       "for the provided (or default) paramNodes:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'warn = FALSE\' in\n",
                       "the control list."))
      }
      reCheck <- setdiff(randomEffectsNodes, reNodesDefault)
      if(length(reCheck)) {
        errorNodes <- paste0(head(reCheck, n = 4), sep = ", ", collapse = ", ")
        if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some randomEffectsNodes provided that look like\n",
                       "they are not needed for Laplace approximation for the\n",
                       "provided (or default) paramNodes:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'warn = FALSE\' in\n",
                       "the control list."))
      }
    }
    if(!reProvided) {
      randomEffectsNodes <- reNodesDefault
    }
    if((!calcProvided) || warn) {
      paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
      reStochDeps <- model$getDependencies(randomEffectsNodes)
      calcNodesDefault <- union(paramStochDeps, reStochDeps)
    }
    if(calcProvided){
      calcNodes <- model$expandNodeNames(calcNodes)
    }
    if(calcProvided  && warn) {
      calcCheck <- setdiff(calcNodesDefault, calcNodes)
      if(length(calcCheck)) {
        errorNodes <- paste0(head(calcCheck, n = 4), sep = ", ", collapse = ", ")
        if(length(calcCheck) > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some model nodes that look like they should be\n",
                       "included in the calcNodes for Laplace approximation over\n",
                       "the provided (or default) randomEffectsNodes:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'warn = FALSE\' in\n",
                       "the control list."))
      }
      calcCheck <- setdiff(calcNodes, calcNodesDefault)
      if(length(calcCheck)){
        errorNodes <- paste0(head(calcCheck, n = 4), sep = ", ", collapse = ", ")
        if(length(calcCheck) > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some calcNodes provided that look like\n",
                       "they are not needed for Laplace approximation over\n",
                       "the provided (or default) randomEffectsNodes:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'warn = FALSE\' in\n",
                       "the control list."))
      }
    }
    if(!calcProvided){
      calcNodes <- calcNodesDefault
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
    ## Create a Laplace nimbleFunctionList
    laplace_nfl <- nimbleFunctionList(Laplace_BASE)
    scalarRENodes <- model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE)
    nre <- length(scalarRENodes)
    ## Record the order of random effects processed internally
    internalRandomEffectsNodes <- NULL
    lenInternalRENodeSets <- NULL
    if(isFALSE(split)) { ## Do all calcNodes in one set
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
      ## Build Laplace
      if(nre > 1) laplace_nfl[[1]] <- nimOneLaplace(model, paramNodes, randomEffectsNodes, calcNodes, innerOptControl, innerOptMethod, innerOptStart)
      else laplace_nfl[[1]] <- nimOneLaplace1D(model, paramNodes, randomEffectsNodes, calcNodes, innerOptControl, "CG", innerOptStart)
    }
    else {## Split calcNodes into conditionally independent sets
      if(isTRUE(split)){
        calcNodesSets <- model$getConditionallyIndependentSets(nodes = calcNodes, givenNodes = paramNodes)
      }
      else if(is.numeric(split)){
        calcNodesSets <- split(calcNodes, split)
      }
      else stop("Invalid value for \'split\' provided in control list")
      num_calcNodesSets <- length(calcNodesSets)
      if(num_calcNodesSets == 0){
        stop("There was a problem determining conditionally independent sets for this model.")
      }
      ## Record calcNodes that do not depend on random effects, e.g. data nodes that depend only on parameter nodes directly
      calcNodesDoNotDepRandomEffects <- NULL
      numLaplace <- 0
      for(i in seq_along(calcNodesSets)){
        ## Need to call getDependencies below because getConditionallyIndependentSets omits deterministic nodes
        these_calcNodes <- model$getDependencies(calcNodesSets[[i]])
        these_reNodes <- intersect(these_calcNodes, randomEffectsNodes)
        ## Check if there are random effects in this set of calcNodes
        if(length(these_reNodes) == 0){
          calcNodesDoNotDepRandomEffects <- c(calcNodesDoNotDepRandomEffects, these_calcNodes)
        }
        else{
          nre_these <- length(model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE))
          internalRandomEffectsNodes <- c(internalRandomEffectsNodes, these_reNodes)
          lenInternalRENodeSets <- c(lenInternalRENodeSets, nre_these)
          numLaplace <- numLaplace + 1
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
          ## If this is the last set of calcNodes, then add all calcNodesDoNotDepRandomEffects, if any, to the calcNodes of the last single Laplace approximation
          if(i == num_calcNodesSets){
            if(length(calcNodesDoNotDepRandomEffects) > 0) these_calcNodes <- c(these_calcNodes, calcNodesDoNotDepRandomEffects)
          }
          ## Build Laplace for each set
          if(nre_these > 1){
            laplace_nfl[[numLaplace]] <- nimOneLaplace(model, paramNodes, these_reNodes, these_calcNodes, innerOptControl, innerOptMethod, innerOptStart)
          }
          else laplace_nfl[[numLaplace]] <- nimOneLaplace1D(model, paramNodes, these_reNodes, these_calcNodes, innerOptControl, "CG", innerOptStart)
        }
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
    paramNodesAsScalars <- model$expandNodeNames(paramNodes, returnScalarComponents = TRUE)
    paramNodesAsScalars_vec <- paramNodesAsScalars
    npar <- length(paramNodesAsScalars)
    if(npar == 1) paramNodesAsScalars_vec <- c(paramNodesAsScalars, "_EXTRA_")
    paramNodesAsScalars_first <- paramNodesAsScalars[1]
    ## setupOutputs(reNodesAsScalars, paramNodesAsScalars)
    ## Automated transformation for parameters
    paramsTransform <- parameterTransform(model, paramNodes)
    pTransform_length <- paramsTransform$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
    ## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    ## Default calculation method for Laplace
    methodID <- 2
    ## Define two nimbleLists for Laplace MLE output
    LaplaceOutputSingleNimbleList <- nimbleList(names = character(1), estimate = double(1), stdError = double(1))
    LaplaceOutputNimbleList <- nimbleList(params = LaplaceOutputSingleNimbleList(),
                                          random = LaplaceOutputSingleNimbleList(),
                                          vcov   = double(2),
                                          scale  = character(0))
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
      one_time_fixes_done <<- TRUE
    },
    Laplace = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      ans <- 0
      for(i in seq_along(laplace_nfl)){
        if(methodID == 1){
          ans <- ans + laplace_nfl[[i]]$Laplace1(p)
        }
        else if(methodID == 2){          
          ans <- ans + laplace_nfl[[i]]$Laplace2(p)
        }
        else if(methodID == 3){
          ans <- ans + laplace_nfl[[i]]$Laplace3(p)
        }
      }
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation w.r.t. parameters
    gr_Laplace = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      ans <- numeric(length = npar)
      for(i in seq_along(laplace_nfl)){
        if(methodID == 1){
          ans <- ans + laplace_nfl[[i]]$gr_Laplace1(p)
        }
        else if(methodID == 2){          
          ans <- ans + laplace_nfl[[i]]$gr_Laplace2(p)
        }
        else if(methodID == 3){
          ans <- ans + laplace_nfl[[i]]$gr_Laplace3(p)
        }
      }
      return(ans)
      returnType(double(1))
    },
    ## Laplace approximation in terms of transformed parameters
    p_transformed_Laplace = function(pTransform = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      p <- paramsTransform$inverseTransform(pTransform)
      ans <- Laplace(p)
      if(is.nan(ans)) ans <- -Inf
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
    derivspInverseTransform = function(pTransform = double(1), order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(pInverseTransform(pTransform), wrt = pTransform_indices, order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    ## Inverse transform random effects to original scale
    reInverseTransform = function(reTrans = double(1)) {
      re <- reTransform$inverseTransform(reTrans)
      return(re)
      returnType(double(1))
    },
    ## Jacobian of the inverse transformation
    derivsreInverseTransform = function(reTrans = double(1), order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(reInverseTransform(reTrans), wrt = reTransform_indices, order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    ## Gradient of the Laplace approximation in terms of transformed parameters
    p_transformed_gr_Laplace = function(pTransform = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      pDerivs <- derivspInverseTransform(pTransform, c(0, 1))
      gr <- gr_Laplace(pDerivs$value) ## pDerivs$value gives original param values
      ans <- (gr %*% pDerivs$jacobian)[1,]
      return(ans)
      returnType(double(1))
    },
    ## Calculate MLEs of parameters
    LaplaceMLE = function(pStart  = double(1, default = Inf),
                          method  = character(0, default = "BFGS"),
                          hessian = logical(0, default = TRUE)) {
      if(any(abs(pStart) == Inf)) pStart <- values(model, paramNodes)
      ## In case parameter nodes are not properly initialized 
      if(any_na(pStart) | any_nan(pStart) | any(abs(pStart)==Inf)) pStartTransform <- rep(0, pTransform_length)
      else pStartTransform <- paramsTransform$transform(pStart)
      ## In case bad start values are provided 
      if(any_na(pStartTransform) | any_nan(pStartTransform) | any(abs(pStartTransform)==Inf)) pStartTransform <- rep(0, pTransform_length)
      optRes <- optim(pStartTransform, p_transformed_Laplace, p_transformed_gr_Laplace,
                      method = method, control = outOptControl, hessian = hessian)
      ## Back transform results to original scale
      optRes$par <- paramsTransform$inverseTransform(optRes$par)
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Optimized random effects given transformed parameter values
    optimRandomEffects = function(pTransform = double(1)){
      p <- pInverseTransform(pTransform)
      raneff <- numeric(nre)
      tmp <- numeric(nre) ## Not sure this is needed. 
      tot <- 0
      for(i in seq_along(laplace_nfl)){
        if(methodID == 1) tmp <- laplace_nfl[[i]]$update_max_inner_logLik_internal(p)
        else tmp <- laplace_nfl[[i]]$update_max_inner_logLik(p)
        numre <- dim(tmp)[1]
        raneff[(tot+1):(tot+numre)] <- tmp
        tot <- tot + numre
      }
      return(raneff)
      returnType(double(1))
    },
    ## Inverse of the negative Hessian of log-likelihood wrt random effects
    inverseNegHess = function(p = double(1), reTransform = double(1)){
      invHess <- matrix(value = 0, nrow = nre, ncol = nre)
      tot <- 0
      for(i in seq_along(laplace_nfl)){
        numre <- lenInternalRENodeSets[i]
        tmp <- laplace_nfl[[i]]$negHess(p, reTransform[(tot+1):(tot+numre)])
        invHess[(tot+1):(tot+numre), (tot+1):(tot+numre)] <- inverse(tmp)
        tot <- tot + numre
      }
      return(invHess)
      returnType(double(2))
    },
    ## Hessian of joint log-likelihood wrt parameters and (transformed) random effects
    hess_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)){
      ans <- matrix(value = 0, nrow = npar, ncol = nre)
      tot <- 0
      for(i in seq_along(laplace_nfl)){
        numre <- lenInternalRENodeSets[i]
        if(methodID == 1){
          tmp <- laplace_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform[(tot+1):(tot+numre)])
        }
        else {
          tmp <- laplace_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re(p, reTransform[(tot+1):(tot+numre)])
        }
        ans[1:npar, (tot+1):(tot+numre)] <- tmp
        tot <- tot + numre
      }
      return(ans)
      returnType(double(2))
    },
    summary = function(LaplaceMLEOutput          = optimResultNimbleList(),
                       originalScale             = logical(0, default = TRUE),
                       calcRandomEffectsStdError = logical(0, default = FALSE), 
                       returnJointCovariance     = logical(0, default = FALSE)){
      if(dim(LaplaceMLEOutput$hessian)[1] == 0) stop("Hessian matrix was not calculated for Laplace MLE")
      ## Output lists
      ans    <- LaplaceOutputNimbleList$new() 
      pres   <- LaplaceOutputSingleNimbleList$new()
      ranres <- LaplaceOutputSingleNimbleList$new()
      ## Parameters
      p <- LaplaceMLEOutput$par
      pTransform <- paramsTransform$transform(p)
      vcov_pTransform <- -inverse(LaplaceMLEOutput$hessian)
      stdErr_pTransform <- sqrt(diag(vcov_pTransform))
      ## Random effects
      optreTransform <- optimRandomEffects(pTransform)
      optre <- reInverseTransform(optreTransform)
      ntot <- npar + nre
      if(returnJointCovariance) {
        ## Inverse of the negative Hessian of log-likelihood wrt transformed random effects at MLEs
        invNegHess <- inverseNegHess(p, optreTransform)
        jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
        jointInvNegHessZero[1:nre, 1:nre] <- invNegHess
        ## Hessian of log-likelihood wrt to params and transformed random effects
        hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
        ## Derivative of inverse transformation for params
        derivspInvTransform  <- derivspInverseTransform(pTransform, c(0, 1))
        JacobpInvTransform   <- derivspInvTransform$jacobian
        ## Jacobian of optimized random effects wrt transformed parameters
        JacobOptreWrtParams <- invNegHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
        jointJacob <- matrix(NA, nrow = ntot, ncol = npar)
        jointJacob[1:nre, 1:npar] <- JacobOptreWrtParams
        jointJacob[(nre+1):ntot, 1:npar] <- diag(npar)
        ## Joint covariance matrix on transformed scale
        vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
        if(originalScale){
          derivsreInvTransform <- derivsreInverseTransform(optreTransform, c(0, 1))
          JacobreInvTransform  <- derivsreInvTransform$jacobian
          JacobJointInvTransform <- matrix(0, nrow = ntot, ncol = ntot)
          JacobJointInvTransform[1:nre, 1:nre] <- JacobreInvTransform
          JacobJointInvTransform[(nre+1):ntot, (nre+1):ntot] <- JacobpInvTransform
          vcov <- JacobJointInvTransform %*% vcov_Transform %*% t(JacobJointInvTransform)
          stdErr_re_p <- sqrt(diag(vcov))
          stdErr_p <- stdErr_re_p[(nre+1):ntot]
          if(calcRandomEffectsStdError){
            ranres$stdError <- stdErr_re_p[1:nre]
          }
          else{
            ranres$stdError <- numeric(0)
          }
          ans$vcov <- vcov
          pres$estimate <- p
          pres$stdError <- stdErr_p
          ranres$estimate <- optre
        }## End of if(originalScale)
        else { ## On transformed scale
          if(calcRandomEffectsStdError){
            stdErr_reTransform <- sqrt(diag(vcov_Transform)[1:nre])
            ranres$stdError <- stdErr_reTransform
          }
          else{
            ranres$stdError <- numeric(0)
          }
          ans$vcov <- vcov_Transform
          pres$estimate <- pTransform
          pres$stdError <- sqrt(diag(vcov_Transform)[(nre+1):ntot])
          ranres$estimate <- optreTransform
        }
      }## End of if(returnJointCovariance)
      else { ## Do not return joint covariance matrix
        ans$vcov <- matrix(nrow = 0, ncol = 0)
        if(originalScale){## On original scale
          pres$estimate <- p
          ranres$estimate <- optre
          if(calcRandomEffectsStdError){
            ## Joint covariance matrix on transform scale
            invNegHess <- inverseNegHess(p, optreTransform)
            jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
            jointInvNegHessZero[1:nre, 1:nre] <- invNegHess
            ## Hessian of log-likelihood wrt to params and transformed random effects
            hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
            ## Derivative of inverse transformation for params
            derivspInvTransform  <- derivspInverseTransform(pTransform, c(0, 1))
            JacobpInvTransform   <- derivspInvTransform$jacobian
            ## Jacobian of optimized random effects wrt transformed parameters
            JacobOptreWrtParams <- invNegHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
            jointJacob <- matrix(NA, nrow = ntot, ncol = npar)
            jointJacob[1:nre, 1:npar] <- JacobOptreWrtParams
            jointJacob[(nre+1):ntot, 1:npar] <- diag(npar)
            ## Join covariance matrix on transformed scale
            vcov_Transform <- jointInvNegHessZero + jointJacob %*% vcov_pTransform %*% t(jointJacob)
            ## Derivatives information
            derivsreInvTransform <- derivsreInverseTransform(optreTransform, c(0, 1))
            JacobreInvTransform  <- derivsreInvTransform$jacobian
            JacobJointInvTransform <- matrix(0, nrow = ntot, ncol = ntot)
            JacobJointInvTransform[1:nre, 1:nre] <- JacobreInvTransform
            JacobJointInvTransform[(nre+1):ntot, (nre+1):ntot] <- JacobpInvTransform
            stdErr <- numeric(ntot)
            for(i in 1:ntot){
              var_i <- (JacobJointInvTransform[i,,drop=FALSE] %*% vcov_Transform %*% t(JacobJointInvTransform[i,,drop=FALSE]))[1,1]
              stdErr[i] <- sqrt(var_i)
            }
            stdErr_p <- stdErr[(nre+1):ntot]
            stdErr_re <- stdErr[1:nre]
            pres$stdError   <- stdErr_p
            ranres$stdError <- stdErr_re
          }## End of if(calcRandomEffectsStdError)
          else { ## Do not calculate standard errors of random effects estimates
            derivspInvTransform  <- derivspInverseTransform(pTransform, c(0, 1))
            JacobpInvTransform   <- derivspInvTransform$jacobian
            stdErr_p <- numeric(npar)
            for(i in 1:npar){
              var_p_i <- (JacobpInvTransform[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobpInvTransform[i,,drop=FALSE]))[1,1]
              stdErr_p[i] <- sqrt(var_p_i)
            }
            pres$stdError <- stdErr_p
            ranres$stdError <- numeric(0)
          }
        }## End of if(originalScale)
        else {## On transformed scale
          pres$estimate <- pTransform
          pres$stdError <- stdErr_pTransform
          ranres$estimate <- optreTransform
          if(calcRandomEffectsStdError){
            invNegHess <- inverseNegHess(p, optreTransform)
            jointInvNegHessZero <- matrix(0, nrow = ntot, ncol = ntot)
            jointInvNegHessZero[1:nre, 1:nre] <- invNegHess
            ## Hessian of log-likelihood wrt to params and transformed random effects
            hessLoglikwrtpre <- hess_logLik_wrt_p_wrt_re(p, optreTransform)
            ## Derivative of inverse transformation for params
            derivspInvTransform  <- derivspInverseTransform(pTransform, c(0, 1))
            JacobpInvTransform   <- derivspInvTransform$jacobian
            ## Jacobian of optimized random effects wrt transformed parameters
            JacobOptreWrtParams <- invNegHess %*% t(hessLoglikwrtpre) %*% JacobpInvTransform
            stdErr_reTransform <- numeric(nre)
            for(i in 1:nre){
              var_reTransform_i <- invNegHess[i, i] + (JacobOptreWrtParams[i,,drop=FALSE] %*% vcov_pTransform %*% t(JacobOptreWrtParams[i,,drop=FALSE]))[1,1]
              stdErr_reTransform[i] <- sqrt(var_reTransform_i)
            }
            ranres$stdError <- stdErr_reTransform
          }
          else{
            ranres$stdError <- numeric(0)
          }
        }
      }
      pres$names <- paramNodesAsScalars_vec
      ranres$names <- reNodesAsScalars_vec
      ans$params <- pres
      ans$random <- ranres
      if(originalScale) ans$scale <- "original"
      else ans$scale <- "transformed"
      return(ans)
      returnType(LaplaceOutputNimbleList())
    }
  ),
  buildDerivs = list(pInverseTransform  = list(),
                     reInverseTransform = list())
)

#' Laplace approximation
#' 
#' Builds a Laplace approximation algorithm for a given NIMBLE model. 
#' 
#' @param model an uncompiled NIMBLE model object.
#' @param paramNodes a character vector of names of parameter nodes in the model; defaults to top-level stochastic nodes.
#' @param randomEffectsNodes a character vector of names of unobserved (latent) nodes to integrate out using the Laplace approximation; defaults to stochastic nodes that depend on \code{paramNodes}.
#' @param calcNodes a character vector of names of nodes for calculating the log-likelihood value; defaults to nodes that depend on \code{randomEffectsNodes} and nodes that only depend on parameter nodes.
#' There may be deterministic nodes between \code{paramNodes} and \code{calcNodes}. These will be included in calculations automatically.
#' @param optimControl a list of control parameters for the inner optimization of Laplace approximation using \code{optim}. Needed only for \code{nimOneLaplace} and \code{nimOneLaplace1D}. See 'Details' of \code{\link{optim}} for further information.
#' @param optimMethod optimization method to be used in \code{optim} for the inner optimization. Needed only for \code{nimOneLaplace} and \code{nimOneLaplace1D}. See 'Details' of \code{\link{optim}}.
#' Currently \code{nimOptim} supports: "\code{Nelder-Mead}", "\code{BFGS}", "\code{CG}", "\code{L-BFGS-B}". By default, method "\code{CG}" is used for \code{nimOneLaplace1D} and "\code{BFGS}" for \code{nimOneLaplace}.
#' @param optimStart choice of start values for the inner optimization. This could be \code{"last"}, \code{"last.best"}, or a vector of user provided values. \code{"last"} means the latest random effects values left in the model will be used. 
#' \code{"last.best"} means the latest random effects values corresponding to currently the largest Laplace likelihood will be used. By default, the initial random effects values will be used for inner optimization.   
#' @param control a named list (for \code{buildLaplace} only) that controls the behavior of the Laplace approximation. See \code{control} section below.
#'
#' @section \code{control} list:
#' 
#' \code{buildLaplace} accepts the following control list elements:
#' \itemize{
#'   \item \code{split}. If TRUE (default), \code{randomEffectsNodes} will be split into conditionally independent sets if possible.
#'         If FALSE, \code{randomEffectsNodes} will be handled as a multivariate block.
#'         If a vector, \code{randomEffectsNodes} will be split by \code{split}(\code{randomEffectsNodes}, \code{control$split}).
#'         The last option allows arbitrary control over how \code{randomEffectsNodes} are blocked.
#'   \item \code{warn}. If TRUE (default), a warning is issued if \code{randomEffectsNodes}/\code{calcNodes} is provided and has extra or missing elements.
#'   \item \code{innerOptimControl}. See \code{optimControl}. 
#'   \item \code{innerOptimMethod}. See \code{optimMethod}.
#'   \item \code{innerOptimStart}. see \code{optimStart}.
#'   \item \code{outOptimControl}. A list of control parameters for maximizing the Laplace log-likelihood using \code{optim}. 
#'         See 'Details' of \code{\link{optim}} for further information.
#' }
#'
#' @section \code{Laplace_BASE}:
#' 
#' Laplace base class, upon which specific Laplace algorithm classes are based by including \code{contains = Laplace_BASE}. This declares a list of nimbleFunctions for a single Laplace approximation.
#' @section \code{nimOneLaplace1D}:
#' 
#' This function is suitable for constructing a single Laplace approximation when \code{randomEffectsNodes} contains only one scalar node.
#' To use this function, one has to accurately provide inputs for all the arguments. 
#' 
#' This function generates an object that comprises a set of methods (functions), each accomplishing one piece of many calculations to obtain the Laplace approximation and its gradient w.r.t. model parameters. 
#' Among these methods, six are most useful to a user:
#' \itemize{
#'   \item \code{Laplace1(p)}. Laplace approximation evaluated at the parameter value \code{p}. This function uses single taping for gradient and Hessian calculations and separate components.
#'   \item \code{Laplace2(p)}. Laplace approximation evaluated at the parameter value \code{p}. This function uses double taping for gradient and Hessian calculations and separate components.
#'   \item \code{Laplace3(p)}. Laplace approximation evaluated at the parameter value \code{p}. This function uses double taping for gradient and Hessian calculations and packs everything together.
#'   \item \code{gr_Laplace1(p)}. Gradient of \code{Laplace1} w.r.t. parameters evaluated at the parameter value \code{p}.
#'   \item \code{gr_Laplace2(p)}. Gradient of \code{Laplace2} w.r.t. parameters evaluated at the parameter value \code{p}.
#'   \item \code{gr_Laplace3(p)}. Gradient of \code{Laplace3} w.r.t. parameters evaluated at the parameter value \code{p}.
#' }
#' 
#' @section \code{nimOneLaplace}:
#' 
#' This function is suitable for constructing a single Laplace approximation when \code{randomEffectsNodes} contains more than one scalar node.
#' To use this function, one has to accurately provide inputs for all the arguments. 
#' 
#' The methods generated by this function are the same as \code{nimOneLaplace1D}. 
#' 
#' @section \code{buildLaplace}:
#' 
#' The main function for constructing the Laplace approximation for a given model. One only needs to provide a NIMBLE model object and then the function
#' will determine inputs for \code{paramNodes}, \code{randomEffectsNodes}, and \code{calcNodes} and then construct the Laplace algorithm. 
#' Default settings are given inside the function and can be changed for all control parameters inside the \code{control} argument. 
#' 
#' The Laplace algorithm object contains a list of functions:
#'\itemize{
#'   \item \code{get_node_name_vec(returnParams)}. Return a list (>1) of names of parameters/random effects nodes, according to \code{returnParams = TRUE/FALSE}.
#'   \item \code{get_node_name_single(returnParams)}. Return the name of a single parameter/random effect node, \code{returnParams = TRUE/FALSE}.
#'   \item \code{set_method(method)}. Set method ID for calculating the Laplace approximation and gradient: 1 (\code{Laplace1}), 2 (\code{Laplace2}, default method), or 3 (\code{Laplace3}).
#'   \item \code{get_method()}. Return the method ID currently used in the algorithm. 
#'   \item \code{one_time_fixes()}. Fix the dimensionality issue if there is only one parameter in the model. This function is called where necessary and users do not need to run this.  
#'   \item \code{Laplace(p)}. Laplace approximation at parameter value \code{p}.
#'   \item \code{gr_Laplace(p)}. Gradient of the Laplace approximation at parameter value \code{p}.
#'   \item \code{p_transformed_Laplace(pTransform)}. Laplace approximation at transformed parameter value \code{pTransform}. 
#'         To make maximizing the Laplace likelihood unconstrained, an automated transformation via \code{\link{parameterTransform}} is performed on any parameters with value constraints.  
#'   \item \code{p_transformed_gr_Laplace(pTransform)}. Gradient of the Laplace approximation (with parameter transformation) w.r.t. transformed parameters, evaluated at transformed parameter value \code{pTransform}.
#'   \item \code{LaplaceMLE(pStart, method, hessian)}. Run maximum likelihood estimation and return results on the transformed scale if any. 
#'         Arguments include \code{pStart}: start value on the original scale; default to parameter values in the model, \code{method}: optimization method used in \code{optim}; default \code{BFGS}, and \code{hessian}: whether calculating the Hessian matrix or not; default to \code{TRUE}.
#'   \item \code{pInverseTransform(pTransform)}. Back transform the transformed parameter value \code{pTransform} to original scale.
#'   \item \code{derivspInverseTransform(pTransform, order)}. Derivative of the inverse transformation w.r.t. transformed parameters at \code{pTransform}. Derivative order is given by \code{order}.
#'   \item \code{reInverseTransform(reTrans)}. Back transform the transformed random effects value \code{reTrans} to original scale.
#'   \item \code{derivsreInverseTransform(reTrans, order)}. Derivative of the inverse transformation w.r.t. transformed random effects at \code{reTrans}. Derivative order is given by \code{order}.
#'   \item \code{optimRandomEffects(pTransform)}. Calculate the optimized random effects given transformed parameter value \code{pTransform}.
#'   \item \code{inverseNegHess(p, reTransform)}. Calculate the inverse of the negative Hessian matrix of the joint log-likelihood w.r.t. transformed random effects, evaluated at parameter value \code{p} and transformed random effects \code{reTransform}.
#'   \item \code{hess_logLik_wrt_p_wrt_re(p, reTransform)}. Calculate the Hessian matrix of the joint log-likelihood w.r.t. parameters and transformed random effects, evaluated at parameter value \code{p} and transformed random effects \code{reTransform}.
#'   \item \code{summary(LaplaceMLEOutput, originalScale, calcRandomEffectsStdError, returnJointCovariance)}. Summarize the maximum likelihood estimation results, given object \code{LaplaceMLEOutput} that is returned by \code{LaplaceMLE}. 
#'        In addition, this function accepts the following arguments:
#'        \itemize{
#'           \item \code{originalScale}. Logical. If TRUE, the function returns results for original parameters and random effects; otherwise, it returns transformed results if transformation is needed. Defaults to TRUE.
#'           \item \code{calcRandomEffectsStdError}. Logical. If TRUE, standard errors of random effects estimators will be calculated. Defaults to FALSE.
#'           \item \code{returnJointCovariance}. Logical. If TRUE, return the joint variance-covariance matrix of the estimators of random effects and parameters. Defaults to FALSE.
#'        } 
#'        This function returns a named list:
#'        \itemize{
#'           \item \code{params}. A list that contains estimates and standard errors of parameters on a specified scale, i.e. original (if \code{originalScale} is TRUE) or transformed.
#'           \item \code{random}. A list that contains estimates of random effects and if required their standard errors, on original/transformed scale. 
#'           \item \code{vcov}. If required (i.e. \code{returnJointCovariance} is TRUE), the joint variance-covariance matrix of the estimators of random effects and parameters, on original/transformed scale. 
#'           \item \code{scale}. \code{original} or \code{transformed}. 
#'        } 
#'}
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
#' # Calculate MLEs on transformed scale
#' MLEres <- CpumpLaplace$LaplaceMLE()
#' # Calculate estimates and standard errors for parameters and random effects on original scale
#' allres <- CpumpLaplace$summary(MLEres, calcRandomEffectsStdError = TRUE)
#' }
#'
NULL
