## NIMBLE AGHQuad approximation
## AGHQuad base class
#' @rdname AGHQuad
#' @export
AGHQuad_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
    AGHQuad1 = function(p = double(1)){
      returnType(double())
    },
    AGHQuad2 = function(p = double(1)){
      returnType(double())
    },
    AGHQuad3 = function(p = double(1)){
      returnType(double())
    },
    gr_AGHQuad1 = function(p = double(1)){
      returnType(double(1))
    },
    gr_AGHQuad2 = function(p = double(1)){
      returnType(double(1))
    },
    gr_AGHQuad3 = function(p = double(1)){
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

## Function to return nodes for Adaptive Gauss-Hermite Quadrature
#' @rdname Rget_AGHQ_nodes
#' @export
Rget_AGHQ_nodes <- function(n = double()) {
	## If n=1 or is miswritten, do Laplace.
	if(n <= 1){	
		z <- c(0,-1)	## For one time fixes
		w_star <- c(sqrt(2*pi), -1)
	}else{
		## Gauss-Hermite Nodes
		gh_nodes <- pracma::gaussHermite(n)
		## Convert to z = sqrt(2)*z, and wgt = exp(x^2)*sqrt(2)
		z <- sqrt(2) * gh_nodes$x
		w_star <- exp(gh_nodes$x^2 + log(gh_nodes$w)) * sqrt(2)
	}
	return(list("z_nodes" = z, "weights" = w_star))
}

## A single AGHQuadrature approximation for only one scalar random effect node
#' @rdname AGHQuad
#' @export
nimOneAGHQuad1D <- nimbleFunction(
  contains = AGHQuad_BASE,
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes, nquad, optimControl, optimMethod, optimStart) {
    ## Check the number of random effects is 1
    nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))
    if(length(nre) != 1) stop("Number of random effects for nimOneAGHQuad1D must be 1.")
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
    
    ## Set up start values for the inner optimization of AGHQuad approximation
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
    
    ## The following are used for caching values and gradient in the AGHQuad3 system
    AGHQuad3_saved_value <- numeric(1)
    AGHQuad3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    AGHQuad3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    ## The following are used for caching values and gradient in the AGHQuad3 system
    max_inner_logLik_saved_par <- as.numeric(c(1, -1))
    max_inner_logLik_saved_value <- numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    ## Record the maximum AGHQuad loglikelihood value for obtaining inner optimization start values
    max_AGHQuad <- -Inf 
    max_AGHQuad_saved_re_value <- as.numeric(c(1, -1))
    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE

	## Adaptive Gauss-Hermite Quadrature Setup:
	AGHQnodes <- Rget_AGHQ_nodes(nquad)
	wgt <- AGHQnodes$weights
	z_node <- AGHQnodes$z_nodes
	middle_node <- -1
	if(nquad %% 2 == 0) middle_node <- (nquad + 1)/2
	
	## Useful cached values for quadrature:
	wi_lik_i <- numeric(nquad + 1)	## Add an extra zero for scalar compilation error.
	sum_wi_lik <- numeric(1)
    quadrature_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    logdetNegHessian <- numeric(1)
	marginal_log_lik <- numeric(1)
	AGHQuad_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
	
    gr_sigmahatwrtre <- numeric(1)
    gr_sigmahatwrtp <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    gr_rehatwrtp <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1)) # double(1)
	gr_QuadSum_value <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
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
      else if(startID == 2) ans <- max_AGHQuad_saved_re_value
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
      max_AGHQuad_saved_re_value <<- fix_one_vec(max_AGHQuad_saved_re_value)
      if(startID == 3) optStart <<- fix_one_vec(optStart)
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        AGHQuad3_saved_gr <<- fix_one_vec(AGHQuad3_saved_gr)
        AGHQuad3_previous_p <<- fix_one_vec(AGHQuad3_previous_p)
        max_inner_logLik_previous_p <<- fix_one_vec(max_inner_logLik_previous_p)
		quadrature_previous_p <<- fix_one_vec(quadrature_previous_p)	## Added by Paul.		
		AGHQuad_saved_gr <<- fix_one_vec(AGHQuad_saved_gr)
		gr_sigmahatwrtp <<- fix_one_vec(gr_sigmahatwrtp)
		gr_rehatwrtp <<- fix_one_vec(gr_rehatwrtp)
		gr_QuadSum_value <<- fix_one_vec(gr_QuadSum_value)
	  }
	  if(nquad == 1) {
	  	wgt <<- fix_one_vec(wgt)
		z_node <<- fix_one_vec(z_node)
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
    ## Solve the inner optimization for AGHQuad approximation
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
        print("Warning: optim does not converge for the inner optimization of AGHQuad approximation")
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
        print("Warning: optim does not converge for the inner optimization of AGHQuad approximation")
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
    ## Put everything (gradient and Hessian) together for AGHQuad3
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
    update_AGHQuad3_with_gr = function(p = double(1), reset = logical(0, default = FALSE)) {
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
      # logdetNegHess_value <- ans$value[ind]
	  logdetNegHessian <<- ans$value[ind]
	  
      ind <- ind + 1
      # chol_negHess <- matrix(ans$value[(ind):(ind + nre*nre - 1)], nrow = nre, ncol = nre)
      negHessValue <- ans$value[ind]
      ind <- ind + 1
      hess_cross_terms <- numeric(value = ans$value[(ind):(ind + npar*1 - 1)], length = npar*1)
      ind <- ind + npar*1
      gr_logdetNegHess_wrt_p_v <- numeric(value = ans$value[(ind):(ind + npar - 1)], length = npar)
      ind <- ind + npar
      gr_logdetNegHess_wrt_re_v <- ans$value[ind]
      
	  #logdetNegHessian <- logdetNegHess(p, reTransform) ## Renamed logdetNegHess_value
	  sigma_hat <- exp(-0.5 * logdetNegHessian)
	
	  ## dre_hat/dp = d^2ll/drep / d^2ll/dre^2
	  gr_rehatwrtp <<- hess_cross_terms/negHessValue

	  ## dsigma_hat/dp (needed at real scale)
	  gr_sigmahatwrtp <<- -0.5*gr_logdetNegHess_wrt_p_v*sigma_hat
	  gr_sigmahatwrtre <<- -0.5*gr_logdetNegHess_wrt_re_v*sigma_hat

	  gr_marg_loglik_wrt_p <- numeric(value = 0, length = length(p))

	  ## AGHQ main call.
	  QuadSum(p) # Caches each log lik evaluations and marginalized log-lik.
	  gr_QuadSum(p, 2)
	  gr_AGHQuad_v <- gr_QuadSum_value  - 0.5 * (gr_logdetNegHess_wrt_p_v + gr_logdetNegHess_wrt_re_v * gr_rehatwrtp)

	  AGHQuad3_saved_gr <<- gr_AGHQuad_v
	  AGHQuad3_saved_value <<- marginal_log_lik
	  return(ans$value)
      returnType(double(1))
    },
    AGHQuad3_update = function(p = double(1)) {
      if(any(p != AGHQuad3_previous_p)) {
        update_AGHQuad3_with_gr(p)
        AGHQuad3_previous_p <<- p
      }
    },
    AGHQuad3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      AGHQuad3_update(p)
      if(AGHQuad3_saved_value > max_AGHQuad) {
        max_AGHQuad <<- AGHQuad3_saved_value
        max_AGHQuad_saved_re_value <<- max_inner_logLik_saved_par
      }
      return(AGHQuad3_saved_value)
      returnType(double())
    },
    gr_AGHQuad3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      AGHQuad3_update(p)
      return(AGHQuad3_saved_gr)
      returnType(double(1))
    },
    ## AGHQuad approximation 2: double tapping with separate components
    AGHQuad2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
	  logdetNegHessian <<- logdetNegHess(p, max_inner_logLik_saved_par)
	  QuadSum(p)
      return(marginal_log_lik)	#Cached and updated by QuadSum
      returnType(double())
    },
    ## AGHQuad approximation 1: single tapping with separate components
    AGHQuad1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
	  logdetNegHessian <<- logdetNegHess(p, max_inner_logLik_saved_par)
	  QuadSum(p)
      return(marginal_log_lik)	#Cached and updated by QuadSum
      returnType(double())
    },	
	## Quadrature sum for both taping methods.
	QuadSum = function(p = double(1)){
	  ## Defaults to using double taping update just in case called somewhere else.
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }

	  ## Track last call of function to match with cached values for gradient.
	  quadrature_previous_p <<- p
	  
	  ## General Quadrature Sum
      reTransform <- max_inner_logLik_saved_par
	  # logdetNegHessian <<- logdetNegHess(p, reTransform)
      maxValue <- max_inner_logLik_saved_value
	  if( nquad == 1 ){	## Laplace Approx.
		ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * 1 * log(2*pi)
		wi_lik_i[1] <<- exp(ans)
	  }else{
        ## AGH approximation:
	    sigma_hat <- exp(- 0.5 * logdetNegHessian)
	    for( i in seq_along(z_node) ){
	      if( middle_node == i ){	## This value we know should be evaluated at the node.
		   wi_lik_i[i] <<- wgt[i]	## exp(maxValue - maxValue) = 1.
		  }else{
	        reTrans_i <- reTransform + sigma_hat*z_node[i]
		    log_lik_i <- joint_logLik(p = p, reTransform = reTrans_i) - maxValue
		    wi_lik_i[i] <<- exp(log_lik_i)*wgt[i]
	      }
	    }
		sum_wi_lik <<- sum(wi_lik_i[1:nquad])	## Needed for gradient.
	    ans <- log(sum_wi_lik) - 0.5 * logdetNegHessian + maxValue
	  }
	  ## Update internally.
   	  marginal_log_lik <<- ans
      if(ans > max_AGHQuad) {
        max_AGHQuad <<- ans
        max_AGHQuad_saved_re_value <<- max_inner_logLik_saved_par
      }
    },
    ## Gradient of the AGHQuad approximation (version 1) w.r.t. parameters
    gr_AGHQuad2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
	  ## Cached transformed RE MLE
	  reTransform <- max_inner_logLik_saved_par

	  ## Need to update cached values for gradient calculation if called 
	  ## before quadrature, which means that we haven't actually calculated
	  if(any(p != quadrature_previous_p ) | marginal_log_lik == 0) {
		logdetNegHessian <<- logdetNegHess(p, reTransform)
		QuadSum(p)
	  }
	  sigma_hat <- exp(-0.5 * logdetNegHessian) ## logdetNegHessian is cached

      ## Double taped gradients.
	  negHessian <- negHess(p, reTransform)[1, 1]
	  grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p(p, reTransform)
	  grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re(p, reTransform)[1]
	  hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re(p, reTransform)[,1]
	  
	  ## dre_hat/dp = d^2ll/drep / d^2ll/dre^2
	  gr_rehatwrtp <<- hesslogLikwrtpre/negHessian
	  ## dsigma_hat/dp (needed at real scale)
	  gr_sigmahatwrtp <<- -0.5*grlogdetNegHesswrtp*sigma_hat
	  gr_sigmahatwrtre <<- -0.5*grlogdetNegHesswrtre*sigma_hat

	  gr_QuadSum(p, 2)
	  AGHQuad_saved_gr <<- gr_QuadSum_value - 0.5 * (grlogdetNegHesswrtp + grlogdetNegHesswrtre * gr_rehatwrtp)
	  return(AGHQuad_saved_gr)
	  returnType(double(1))
	},
    ## Gradient of the AGHQuad approximation (version 1) w.r.t. parameters
    gr_AGHQuad1 = function(p = double(1)){
		if(!one_time_fixes_done) one_time_fixes()
		if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
			update_max_inner_logLik_internal(p)
		}
	  ## Cached transformed RE MLE
	  reTransform <- max_inner_logLik_saved_par

	  ## Need to update cached values for gradient calculation if called 
	  ## before quadrature, which means that we haven't actually calculated
	  if(any(p != quadrature_previous_p ) | marginal_log_lik == 0) {
		logdetNegHessian <<- logdetNegHess(p, reTransform)
		QuadSum(p)
	  }
	  sigma_hat <- exp(-0.5 * logdetNegHessian) ## logdetNegHessian is cached

      ## Double taped gradients.
	  negHessian <- negHess_internal(p, reTransform)[1, 1]
	  grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p_internal(p, reTransform)
	  grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re_internal(p, reTransform)[1]
	  hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform)[,1]
	  
	  ## dre_hat/dp = d^2ll/drep / d^2ll/dre^2
	  gr_rehatwrtp <<- hesslogLikwrtpre/negHessian
	  ## dsigma_hat/dp (needed at real scale)
	  gr_sigmahatwrtp <<- -0.5*grlogdetNegHesswrtp*sigma_hat
	  gr_sigmahatwrtre <<- -0.5*grlogdetNegHesswrtre*sigma_hat

	  gr_QuadSum(p, 1)
	  AGHQuad_saved_gr <<- gr_QuadSum_value - 0.5 * (grlogdetNegHesswrtp + grlogdetNegHesswrtre * gr_rehatwrtp)
	  return(AGHQuad_saved_gr)
	  returnType(double(1))
	},
	## General Quadrature Gradient Summation for all 3 methods.
	## Must be called internally as it does not safe check for previous QuadSum call etc.
	gr_QuadSum = function(p = double(1), method = double())
	{
	  gr_lik_wrt_p <- numeric(value = 0, length = length(p))

      reTransform <- max_inner_logLik_saved_par
	  sigma_hat <- exp(-0.5 * logdetNegHessian)
	  
	  if( method == 2 ) {
		gr_jointlogLikwrtp <- gr_joint_logLik_wrt_re(p, reTransform)
	  }else {
		gr_jointlogLikwrtp <- gr_joint_logLik_wrt_re_internal(p, reTransform)
	  }
	  
	  if( nquad == 1 ){	## Laplace Approx from existing code
	    ## dll/dp + log(sigma)/dp + d^2ll/dre^2 * dre/dp
		gr_QuadSum_value <<- gr_jointlogLikwrtp
	  }else{
	    for( i in seq_along(z_node) ){
	      if( middle_node == i ){	## Try and avoid this set up. Ask Chris where I should put this then?
		    gr_lik_wrt_p <- gr_lik_wrt_p + wi_lik_i[i]*gr_jointlogLikwrtp
		  }else {			
	        reTrans_i <- reTransform + sigma_hat*z_node[i]
		    ## Chain Rule: dll/dre * ( dre_hat/dp + dsigma_hat/dp*z_i )##Improve name
		    ## dll/dp
		    ## Double Taping Version
		    if( method == 2 ) {
		      gr_logLikwrtrewrtp_i <- gr_joint_logLik_wrt_re(p, reTrans_i)[1] * 
				( (1 + gr_sigmahatwrtre*z_node[i]) * gr_rehatwrtp  +  gr_sigmahatwrtp*z_node[i] )
		      gr_logLikewrtp_i <- gr_joint_logLik_wrt_p(p, reTrans_i)	## single call for both RE and p.					
		    }else {
		      ## Single Taping
		      gr_logLikwrtrewrtp_i <- gr_joint_logLik_wrt_re_internal(p, reTrans_i)[1] * 
				( (1 + gr_sigmahatwrtre*z_node[i]) * gr_rehatwrtp  +  gr_sigmahatwrtp*z_node[i] )
		      gr_logLikewrtp_i <- gr_joint_logLik_wrt_p_internal(p, reTrans_i)	## single call for both RE and p.										
		    }
		    ## The weighted gradient for the ith sum.
		    gr_lik_wrt_p <- gr_lik_wrt_p + wi_lik_i[i]*(gr_logLikewrtp_i +  gr_logLikwrtrewrtp_i)
	      }
	    }	
	    gr_QuadSum_value <<- gr_lik_wrt_p/sum_wi_lik
 	  }
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
) ## End of nimOneAGHQuad1D

#' Process model to organize nodes for AGHQuad approximation
#' @export
setupAGHQuadNodes <- function(model, paramNodes, randomEffectsNodes, calcNodes,
                              calcNodesNoAGHQuad, allowNonPriors = FALSE, split = TRUE,
                              warn = TRUE) {
  paramProvided         <- !missing(paramNodes)
  reProvided            <- !missing(randomEffectsNodes)
  calcProvided          <- !missing(calcNodes)
  calcNoAGHQuadProvided <- !missing(calcNodesNoAGHQuad)
  # We may need to use determ and stochastic dependencies of parameters multiple times below
  # Define these to avoid repeated computation
  paramDetermDeps <- character(0)
  paramStochDeps  <- character(0)
  paramDetermDepsCalculated <- FALSE
  paramStochDepsCalculated  <- FALSE
  
  # 1. Default parameters are stochastic top-level nodes
  #    unless allowNonPriors is TRUE, in which case they are all top-level nodes
  if(!paramProvided) {
    allTopNodes <- model$getNodeNames(topOnly = TRUE, includeRHSonly = TRUE)
    stochTopNodes <- allTopNodes[which(model$getNodeType(allTopNodes)=="stoch")]
    if(!allowNonPriors) { 
      paramNodes <- stochTopNodes
    }
    else { 
      allRHSnodes <- allTopNodes[which(model$getNodeType(allTopNodes)=="RHSonly")]
      ## Select all RHSonly nodes with stochastic dependencies
      allRHSwithStochDeps <- character(0)
      numAllRHSnodes <- length(allRHSnodes)
      if(numAllRHSnodes){
        for(i in 1:numAllRHSnodes) {
          if(length(model$getDependencies(allRHSnodes[i], stochOnly = TRUE))) 
            allRHSwithStochDeps <- c(allRHSwithStochDeps, allRHSnodes[i])
        }
      }
      paramNodes <- c(stochTopNodes, allRHSwithStochDeps)
    }
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
                     "like they should be included in randomEffectsNodes for AGHQuad approximation\n",
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
                     "they are not needed for AGHQuad approximation for the\n",
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
  #    Note that the internal AGHQuad nimbleFunctions (one for each conditionally independent set)
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
                     "included in the calcNodes for AGHQuad approximation because\n",
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
    # Check extra determ nodes
    if(numDetermCalcCheck){
      paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      paramDetermDepsCalculated <- TRUE
      for(i in 1:numDetermCalcCheck){
        if(!(determCalcCheck[i] %in% paramDetermDeps) || 
           !(any(m$getDependencies(determCalcCheck[i], self = FALSE) %in% calcNodesDefault))){
          errorNodes <- c(errorNodes, determCalcCheck[i])
        }
      }
    }
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      warning(paste0("There are some calcNodes provided that look like\n",
                     "they are not needed for AGHQuad approximation over\n",
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
  #    that are not involved in the marginalization done by AGHQuad.
  #    Default is a bit complicated: All dependencies from paramNodes to
  #    stochastic nodes that are not part of calcNodes. Note that calcNodes
  #    does not necessarily contain deterministic nodes between paramNodes and
  #    randomEffectsNodes. We don't want to include those in calcNodesNoAGHQuad.
  #    (A deterministic that is needed for both calcNodes and calcNodesNoAGHQuad should be included.)
  #    So we have to first do a setdiff on stochastic nodes and then fill in the
  #    deterministics that are needed.
  if(!calcNoAGHQuadProvided || warn) { 
    if(!paramStochDepsCalculated){
      paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
      paramStochDepsCalculated <- TRUE
    }
    calcNodesNoAGHQuadDefault <- setdiff(paramStochDeps, calcNodes)
  }
  if(calcNoAGHQuadProvided){
    if(is.null(calcNodesNoAGHQuad) || isFALSE(calcNodesNoAGHQuad)) calcNodesNoAGHQuad <- character(0)
    calcNodesNoAGHQuad <- model$expandNodeNames(calcNodesNoAGHQuad, sort = TRUE)
    if((length(calcNodesNoAGHQuad) > 0) && !any(model$getNodeType(calcNodesNoAGHQuad)=="stoch")){
      warning("There are no stochastic nodes in the calcNodesNoAGHQuad provided for AGHQuad.")
    }
  }
  if(!calcNoAGHQuadProvided){
    calcNodesNoAGHQuad <- calcNodesNoAGHQuadDefault
  }
  if(calcNoAGHQuadProvided && warn) {
    calcNoAGHQuadCheck <- setdiff(calcNodesNoAGHQuadDefault, calcNodesNoAGHQuad)
    if(length(calcNoAGHQuadCheck)) {
      # We only check missing stochastic nodes; determ nodes will be added below
      missingStochNodesInds <- which((model$getNodeType(calcNoAGHQuadCheck)) == "stoch")
      numMissingStochNodes <- length(missingStochNodesInds)
      missingStochNodes <- calcNoAGHQuadCheck[missingStochNodesInds]
      if(numMissingStochNodes){
        errorNodes <- paste0(head(missingStochNodes, n = 4), sep = "", collapse = ", ")
        if(numMissingStochNodes > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some model nodes (stochastic) that look like they should be\n",
                       "included in the calcNodesNoAGHQuad for exact likelihood calculation:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'warn = FALSE\' in the control list."))
      }
    }
    # Check redundant stochastic nodes
    calcNoAGHQuadCheck <- setdiff(calcNodesNoAGHQuad, calcNodesNoAGHQuadDefault)
    stochCalcNoAGHQuadCheck <- calcNoAGHQuadCheck[model$getNodeType(calcNoAGHQuadCheck)=="stoch"]
    # Check redundant determ nodes
    determCalcNoAGHQuadCheck <- setdiff(calcNoAGHQuadCheck, stochCalcNoAGHQuadCheck)
    numDetermCalcNoAGHQuadCheck <- length(determCalcNoAGHQuadCheck)
    errorNodes <- character(0)
    if(numDetermCalcNoAGHQuadCheck){
      if(!paramDetermDepsCalculated) {
        paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
        paramDetermDepsCalculated <- TRUE
      }
      for(i in 1:numDetermCalcNoAGHQuadCheck){
        if(!(determCalcNoAGHQuadCheck[i] %in% paramDetermDeps) || 
           !(any(m$getDependencies(determCalcNoAGHQuadCheck[i], self = FALSE) %in% calcNodesNoAGHQuadDefault))){
          errorNodes <- c(errorNodes, determCalcNoAGHQuadCheck[i])
        }
      }
    }
    errorNodes <- c(stochCalcNoAGHQuadCheck, errorNodes)
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      warning(paste0("There are some nodes provided in calcNodesNoAGHQuad that look like\n",
                     "they are not needed for exact likelihood calculation:\n",
                     outErrorNodes, "\n",
                     "To silence this warning, include \'warn = FALSE\' in the control list."))
    }
  }
  # Check and add necessary (upstream) deterministic nodes into calcNodesNoAGHQuad
  # This ensures that deterministic nodes between paramNodes and calcNodesNoAGHQuad are used.
  num_calcNodesNoAGHQuad <- length(calcNodesNoAGHQuad)
  if(num_calcNodesNoAGHQuad > 0){
    if(!paramDetermDepsCalculated) {
      paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE)
      paramDetermDepsCalculated <- TRUE
    }
    numParamDetermDeps <- length(paramDetermDeps)
    if(numParamDetermDeps > 0) {
      keep_paramDetermDeps <- logical(numParamDetermDeps)
      for(i in seq_along(paramDetermDeps)) {
        nextDeps <- model$getDependencies(paramDetermDeps[i])
        keep_paramDetermDeps[i] <- any(nextDeps %in% calcNodesNoAGHQuad)
      }
      paramDetermDeps <- paramDetermDeps[keep_paramDetermDeps]
    }
    calcNodesNoAGHQuad <- model$expandNodeNames(c(paramDetermDeps, calcNodesNoAGHQuad), sort = TRUE)
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
       calcNodesNoAGHQuad = calcNodesNoAGHQuad,
       givenNodes = givenNodes,
       randomEffectsSets = reSets
       )
}

## Main function for AGHQuad approximation
#' @rdname AGHQuad 
#' @export
buildAGHQuad <- nimbleFunction(
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes, calcNodesNoAGHQuad, 
                   nquad, control = list()) {
    if(is.null(control$split)) split <- TRUE else split <- control$split
    if(is.null(control$warn))   warn <- TRUE else  warn <- control$warn
    if(is.null(control$allowNonPriors)) allowNonPriors <- FALSE else  allowNonPriors <- control$allowNonPriors

	nquadProvided  <- !missing(nquad)    
    if(!nquadProvided) {
	  nquad <- 5					  
      nimCat('  [Note] Assuming 5 quadrature nodes.\n')
	}
    AGHQuadNodes <- NULL
    if(!missing(paramNodes)) {
      if(is.list(paramNodes)) {
        # The user called setupAGHQuadNodes and provided a list of that format to paramNodes.
        AGHQuadNodes <- paramNodes
      }
    }
    if(is.null(AGHQuadNodes)) {
      AGHQuadNodes <- setupAGHQuadNodes(model = model, paramNodes = paramNodes,
                                        randomEffectsNodes = randomEffectsNodes,
                                        calcNodes = calcNodes,
                                        calcNodesNoAGHQuad = calcNodesNoAGHQuad,
                                        allowNonPriors = allowNonPriors,
                                        split = split,
                                        warn = warn)
    }
    paramNodes <- AGHQuadNodes$paramNodes
    randomEffectsNodes <- AGHQuadNodes$randomEffectsNodes
    calcNodes <- AGHQuadNodes$calcNodes
    calcNodesNoAGHQuad <- AGHQuadNodes$calcNodesNoAGHQuad
    num_calcNodesNoAGHQuad <- length(calcNodesNoAGHQuad)
    # AGHQuadNodes$randomEffectsSets will be extracted below if needed

    if(length(calcNodesNoAGHQuad)) {
      exactLoglik_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = paramNodes, calcNodes = calcNodesNoAGHQuad)
      exactLoglik_updateNodes   <- exactLoglik_derivsInfo$updateNodes
      exactLoglik_constantNodes <- exactLoglik_derivsInfo$constantNodes
    }
    else { ## calcNodesNoAGHQuad is empty
      exactLoglik_updateNodes   <- character(0)
      exactLoglik_constantNodes <- character(0)
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
    
    ## Create a AGHQuad nimbleFunctionList
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
        if(nre > 1) {
			nimCat('  [Note] Adaptive Gauss-Hermite quadrature is not implemented for multivariate integration.\n Defaulting to Laplace Approximation.')
			AGHQuad_nfl[[1]] <- nimOneLaplace(model, paramNodes, randomEffectsNodes, calcNodes, innerOptControl, innerOptMethod, innerOptStart)
		}else AGHQuad_nfl[[1]] <- nimOneAGHQuad1D(model, paramNodes, randomEffectsNodes, calcNodes, nquad, innerOptControl, "CG", innerOptStart)
      }
      else {## Split randomEffectsNodes into conditionally independent sets
        reSets <- AGHQuadNodes$randomEffectsSets
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
		 	nimCat('  [Note] Adaptive Gauss-Hermite quadrature is not implemented for multivariate integration.\n Defaulting to Laplace Approximation.')
            AGHQuad_nfl[[i]] <- nimOneLaplace(model, paramNodes, these_reNodes, these_calcNodes, innerOptControl, innerOptMethod, innerOptStart)
          }
          else AGHQuad_nfl[[i]] <- nimOneAGHQuad1D(model, paramNodes, these_reNodes, these_calcNodes, nquad, innerOptControl, "CG", innerOptStart)
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
      if(num_calcNodesNoAGHQuad == 0) 
        stop("Both calcNodesNoAGHQuad and randomEffectsNodes are empty for AGHQuad for the given model.")
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
    ## Define two nimbleLists for AGHQuad MLE output
    AGHQuadOutputSingleNimbleList <- nimbleList(names = character(1), estimate = double(1), stdError = double(1))
    AGHQuadOutputNimbleList <- nimbleList(params = AGHQuadOutputSingleNimbleList(),
                                          random = AGHQuadOutputSingleNimbleList(),
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
      if(nre == 0) print("AGHQuad approximation is not needed for the given model: no random effects")
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
    ## Exact log-likelihood in terms of original parameters
    logLikNoAGHQuad = function(p = double(1)) {
      if(num_calcNodesNoAGHQuad == 0) stop("calcNodesNoAGHQuad is empty: there is no exact likelihood component for the model")
      values(model, paramNodes) <<- p
      ans <- model$calculate(calcNodesNoAGHQuad) 
      return(ans)
      returnType(double())
    },
    ## Gradient of the exact log-likelihood w.r.t parameters
    gr_logLikNoAGHQuad_internal = function(p = double(1)) {
      if(num_calcNodesNoAGHQuad == 0) stop("calcNodesNoAGHQuad is empty: there is no exact likelihood component for the model")
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(logLikNoAGHQuad(p), wrt = p_indices, order = 1, model = model, 
                    updateNodes = exactLoglik_updateNodes, constantNodes = exactLoglik_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    gr_logLikNoAGHQuad = function(p = double(1)) {
      if(num_calcNodesNoAGHQuad == 0) stop("calcNodesNoAGHQuad is empty: there is no exact likelihood component for the model")
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(gr_logLikNoAGHQuad_internal(p), wrt = p_indices, order = 0, model = model,
                    updateNodes = exactLoglik_updateNodes, constantNodes = exactLoglik_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
	Laplace = function(p = double(1)) {
		if(length(p) != npar) stop("  [Warning] Calculation requires values for all parameters.")
		if(nquad != 1) {
			print('  [Note] Applying Adaptive Gauss-Hermite Quadrature with ', nquad, ' nodes.')
		}
		return(calcMargLogLik(p))
		returnType(double())
	},
    ## AGHQuad approximation in terms of original parameters
    calcMargLogLik = function(p = double(1)) {
      if(length(p) != npar) stop("  [Warning] Calculation requires values for all parameters.")
      if(!one_time_fixes_done) one_time_fixes()
      if(num_calcNodesNoAGHQuad > 0) ans <- logLikNoAGHQuad(p)
      else ans <- 0
      if(nre > 0){
        for(i in seq_along(AGHQuad_nfl)){
          if(methodID == 1) ans <- ans + AGHQuad_nfl[[i]]$AGHQuad1(p)
          else if(methodID == 2) ans <- ans + AGHQuad_nfl[[i]]$AGHQuad2(p)
          else ans <- ans + AGHQuad_nfl[[i]]$AGHQuad3(p)
        }
      }
      return(ans)
      returnType(double())
    },
    ## Gradient of the AGHQuad approximation w.r.t. parameters
    gr_margLogLik = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      if(num_calcNodesNoAGHQuad > 0) ans <- gr_logLikNoAGHQuad(p)
      else ans <- numeric(length = npar)
      if(nre > 0){
        for(i in seq_along(AGHQuad_nfl)){
          if(methodID == 1) ans <- ans + AGHQuad_nfl[[i]]$gr_AGHQuad1(p)
          else if(methodID == 2) ans <- ans + AGHQuad_nfl[[i]]$gr_AGHQuad2(p)
          else ans <- ans + AGHQuad_nfl[[i]]$gr_AGHQuad3(p)
        }
      }
      return(ans)
      returnType(double(1))
    },
    ## AGHQuad approximation in terms of transformed parameters
    p_transformed_margLogLik = function(pTransform = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      p <- paramsTransform$inverseTransform(pTransform)
      ans <- calcMargLogLik(p) 
      if(is.nan(ans) | is.na(ans)) ans <- -Inf
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
      if(nre == 0) stop("No random effects in the model")
      re <- reTransform$inverseTransform(reTrans)
      return(re)
      returnType(double(1))
    },
    ## Jacobian of the inverse transformation
    derivsreInverseTransform = function(reTrans = double(1), order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      if(nre == 0) stop("No random effects in the model")
      ans <- derivs(reInverseTransform(reTrans), wrt = reTransform_indices, order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    ## Gradient of the AGHQuad approximation in terms of transformed parameters
    p_transformed_gr_margLogLik = function(pTransform = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      pDerivs <- derivspInverseTransform(pTransform, c(0, 1))
      gr <- gr_margLogLik(pDerivs$value) ## pDerivs$value gives original param values
      ans <- (gr %*% pDerivs$jacobian)[1,] 
      return(ans)
      returnType(double(1))
    },
	## Pointer for old function name.
	LaplaceMLE = function(pStart  = double(1, default = Inf),
                          method  = character(0, default = "BFGS"),
                          hessian = logical(0, default = TRUE)) {
		if(nquad != 1) {
			print('  [Note] Applying Adaptive Gauss-Hermite Quadrature with ', nquad, ' nodes.')
		}
		return(calculateMLE(pStart, method, hessian))
        returnType(optimResultNimbleList())
	},
    ## Calculate MLEs of parameters
    calculateMLE = function(pStart  = double(1, default = Inf),
                          method  = character(0, default = "BFGS"),
                          hessian = logical(0, default = TRUE)) {
      if(any(abs(pStart) == Inf)) pStart <- values(model, paramNodes)
      ## In case parameter nodes are not properly initialized 
      if(any_na(pStart) | any_nan(pStart) | any(abs(pStart)==Inf)) pStartTransform <- rep(0, pTransform_length)
      else pStartTransform <- paramsTransform$transform(pStart)
      ## In case bad start values are provided 
      if(any_na(pStartTransform) | any_nan(pStartTransform) | any(abs(pStartTransform)==Inf)) pStartTransform <- rep(0, pTransform_length)
      optRes <- optim(pStartTransform, p_transformed_margLogLik, p_transformed_gr_margLogLik, method = method, control = outOptControl, hessian = hessian)
      if(optRes$convergence != 0) 
        print("Warning: optim has a non-zero convergence code: ", optRes$convergence, ".\nThe control parameters of optim can be adjusted in the control argument of buildAGHQuad via list(outOptimControl = list()).")
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
    inverseNegHess = function(p = double(1), reTransform = double(1)){
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
    summary = function(AGHQuadMLEOutput          = optimResultNimbleList(),
                       originalScale             = logical(0, default = TRUE),
                       calcRandomEffectsStdError = logical(0, default = FALSE), 
                       returnJointCovariance     = logical(0, default = FALSE)){
      if(dim(AGHQuadMLEOutput$hessian)[1] == 0) stop("Hessian matrix was not calculated for AGHQuad MLE")
      ## Output lists
      ans <- AGHQuadOutputNimbleList$new() 
      pres <- AGHQuadOutputSingleNimbleList$new()
      ranres <- AGHQuadOutputSingleNimbleList$new()
      ## Parameters
      p <- AGHQuadMLEOutput$par
      pTransform <- paramsTransform$transform(p)
      vcov_pTransform <- -inverse(AGHQuadMLEOutput$hessian)
      stdErr_pTransform <- sqrt(diag(vcov_pTransform))
      if(nre == 0) { ## No random effects
        ranres$estimate <- numeric(0)
        ranres$stdError <- numeric(0)
        if(originalScale){
          derivspInvTransform  <- derivspInverseTransform(pTransform, c(0, 1))
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
          pres$estimate <- p
          pres$stdError <- stdErr_p
        }
        else {
          pres$estimate <- pTransform
          pres$stdError <- stdErr_pTransform
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
      }
      pres$names <- paramNodesAsScalars_vec
      ranres$names <- reNodesAsScalars_vec
      ans$params <- pres
      ans$random <- ranres
      if(originalScale) ans$scale <- "original"
      else ans$scale <- "transformed"
      return(ans)
      returnType(AGHQuadOutputNimbleList())
    }
  ),
  buildDerivs = list(pInverseTransform  = list(),
                     reInverseTransform = list(),
                     logLikNoAGHQuad = list(),
                     gr_logLikNoAGHQuad_internal = list()
                     )
)

#' AGHQuad approximation
#' 
#' Builds a AGHQuad approximation algorithm for a given NIMBLE model. 
#' 
#' @param model a NIMBLE model object, such as returned by \code{nimbleModel}.
#'   The model must have automatic derivatives (AD) turned on, using
#'   \code{buildDerivs=TRUE} in \code{nimbleModel}.
#' @param paramNodes a character vector of names of parameter nodes in the
#'   model; defaults to top-level stochastic nodes, as determined by
#'   \code{model$getNodeNames(topOnly=TRUE, stochOnly=TRUE)}. If
#'   \code{allowNonPriors} is \code{TRUE}, top-level determinisic nodes are also
#'   treated as parameters (i.e. \code{stochOnly=FALSE}). Alternatively,
#'   \code{paramNodes} can be a list in the format returned by
#'   \code{setupAGHQuadNodes}, in which case \code{randomEffectsNodes},
#'   \code{calcNodes}, and \code{calcNodesNoAGHQuad} are not needed (and will be
#'   ignored).
#' @param randomEffectsNodes a character vector of names of unobserved (latent)
#'   nodes to marginalize (integrate) over using the AGHQuad approximation;
#'   defaults to latent stochastic nodes that depend on \code{paramNodes}.
#' @param calcNodes a character vector of names of nodes for calculating the
#'   integrand for AGHQuad approximation; defaults to nodes that depend on
#'   \code{randomEffectsNodes}, as determined by
#'   \code{model$geteDependencies(randomEffectsNodes)} (which will include
#'   \code{randomEffectsNodes}). There may be deterministic nodes between
#'   \code{paramNodes} and \code{randomEffectsNodes}. These will be included in
#'   calculations automatically and thus do not need to be included in
#'   \code{calcNodes} (but there is no problem if they are).
#' @param calcNodesNoAGHQuad a character vector of names of nodes for
#'   calculating terms in the log-likelihood that do not depend on any
#'   \code{randomEffectsNodes}, and thus are not part of the marginalization,
#'   but should be included for purposes of finding the MLE. This defaults to
#'   stochastic nodes that depend on \code{paramNodes} but are not part of and
#'   do not depend on \code{randomEffectsNodes}. There may be deterministic
#'   nodes between \code{paramNodes} and \code{calcNodesNoAGHQuad}. These will
#'   be included in calculations automatically and thus do not need to be
#'   included in \code{calcNodesNoAGHQuad} (but there is no problem if they
#'   are).
#' @param optimControl a list of control parameters for the inner optimization
#'   (of randomEffectsNodes) of AGHQuad approximation using \code{optim}. This
#'   is used in the internal nimbleFunctions \code{nimOneAGHQuad} and
#'   \code{nimOneAGHQuad1D}. See 'Details' of \code{\link{optim}} for further
#'   information.
#' @param optimMethod optimization method to be used in \code{optim} for the
#'   inner optimization. This is used in the internal nimbleFunctions
#'   \code{nimOneAGHQuad} and \code{nimOneAGHQuad1D}. See 'Details' of
#'   \code{\link{optim}}. Currently \code{nimOptim} supports:
#'   "\code{Nelder-Mead}", "\code{BFGS}", "\code{CG}", "\code{L-BFGS-B}". By
#'   default, method "\code{CG}" is used for \code{nimOneAGHQuad1D} and
#'   "\code{BFGS}" for \code{nimOneAGHQuad}.
#' @param optimStart choice of starting values for the inner optimization. This
#'   could be \code{"last"}, \code{"last.best"}, or a vector of user provided
#'   values. \code{"last"} means the most recent random effects values left in
#'   the model will be used. When finding the MLE, the most recent values will
#'   be the result of the most recent inner optimization, from the previous
#'   parameter values. \code{"last.best"} means the random effects values
#'   corresponding to the largest AGHQuad likelihood (from any call to the
#'   `AGHQuad` method, including during an MLE search) will be used (even if it
#'   was not the most recent AGHQuad likelihood). By default, the initial random
#'   effects values will be used for inner optimization.
#' @param control a named list (for \code{buildAGHQuad} only) that controls the
#'   behavior of the AGHQuad approximation. See \code{control} section below.
#'
#' @section \code{buildAGHQuad}:
#'
#' This is the main function for constructing the AGHQuad approximation for a
#' given model. One only needs to provide a NIMBLE model object and then the
#' function will construct the pieces necessary for AGHQuad approximation to
#' marginalize over all latent states (aka random effects) in a model. To do so,
#' it will determine default values for \code{paramNodes},
#' \code{randomEffectsNodes}, \code{calcNodes}, and \code{calcNodesNoAGHQuad} as
#' described above.
#'
#' The default values are obtained by calling \code{setupAGHQuadNodes}, whose
#' arguments match those here (except for a few arguments which are taken from
#' control list elements here). One can call that function to see exactly how
#' nodes will be arranged for AGHQuad approximation. One can also call it,
#' customize the returned list, and then provide that to \code{buildAGHQuad} as
#' \code{paramNodes}.
#'
#' If any \code{paramNodes} (parameters) or \code{randomEffectsNodes} (random
#' effects / latent states) have constraints on the range of valid values
#' (because of the distribution they follow), they will be used on a transformed
#' scale determined by \code{parameterTransform}. This means the AGHQuad
#' approximation itself will be done on the transformed scale for random effects
#' and finding the MLE will be done on the transformed scale for parameters. For
#' parameters, any prior distributions are not included in calculations, but
#' they are used to determine valid parameter ranges. For example, if
#' \code{sigma} is a standard deviation, declare it with a prior such as
#' \code{sigma ~ dhalfflat()} to indicate that it must be greater than 0.
#'
#' The object returned by \code{buildAGHQuad} is a nimbleFunction object with
#' numerous methods (functions). The most useful ones are:
#'
#' \itemize{
#'
#' \item \code{AGHQuad(p)}. AGHQuad approximation to the marginal log-likelihood
#'       function at parameter value \code{p}, which should match the order of
#'       \cope{paramNodes}. For any non-scalar nodes in \code{paramNodes}, the
#'       order within the node is column-major (which can be seen for R objects
#'       using \code{as.numeric}).
#'
#' \item \code{AGHQuadMLE(pStart, method, hessian)}. Find the maximum likelihood
#'         estimates of the AGHQuad-approximated marginal likelihood. Arguments
#'         include \code{pStart}: initial parameter values (defaults to
#'         parameter values currently in the model); \code{method}: (outer)
#'         optimization method to use in \code{optim} (defaults to "BFGS"); and
#'         \code{hessian}: whether to calculate and return the Hessian matrix
#'         (defaults to \code{TRUE}). Second derivatives in the Hessian are
#'         determined by finite differences of the gradients obtained by
#'         automatic differentiation (AD).
#'
#' \item \code{summary(AGHQuadMLEOutput, originalScale,
#'        calcRandomEffectsStdError, returnJointCovariance)}. Summarize the
#'        maximum likelihood estimation results, given object
#'        \code{AGHQuadMLEOutput} that was returned by \code{AGHQuadMLE}. The
#'        summary can include a covariance matrix for the parameters, the random
#'        effects, or both, and these can be returned on the original parameter
#'        scale or on the potentially transformed scale(s) used in estimation.
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
#'        \code{summary} function returns a named list, with elements:
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
#'           \item \code{scale}. \code{original} or \code{transformed}, the
#'        scale on which results were requested.
#'
#'        }
#'
#'     }
#'
#' Additional methods to access or control more details of the AGHQuad approximation include:
#'
#' \itemize{
#'
#'   \item \code{get_node_name_vec(returnParams)}. Return a vector (>1) of names
#'   of parameters/random effects nodes, according to \code{returnParams =
#'   TRUE/FALSE}. Use this is there is more than one parameter.
#'
#'   \item \code{get_node_name_single(returnParams)}. Return the name of a
#'   single parameter/random effect node, \code{returnParams = TRUE/FALSE}. Use
#'   this if there is only one parameter.
#'
#'   \item \code{set_method(method)}. Set method ID for calculating the AGHQuad
#'   approximation and gradient: 1 (\code{AGHQuad1}), 2 (\code{AGHQuad2},
#'   default method), or 3 (\code{AGHQuad3}). See below for more details. Users
#'   wanting to explore efficiency can try switching from method 2 (default) to
#'   methods 1 or 3 and comparing performance. The first AGHQuad approximation
#'   with each method will be (much) slower than subsequent AGHQuad
#'   approximations.
#'
#'   \item \code{get_method()}. Return the current method ID.
#'
#'   \item \code{one_time_fixes()}. Users do not need to run this. Is is called
#'   when necessary internally to fix dimensionality issues if there is only
#'   one parameter in the model.
#'
#'   \item \code{logLikNoAGHQuad(p)}. Calculate the \code{calcNodesNoAGHQuad}
#'   nodes, which returns the log-likelihood of the parts of the model that are
#'   not included in the AGHQuad approximation. \code{p} is the vector of
#'   parameter values in the order of \code{paramNames}.
#'
#'   \item \code{gr_logLikNoAGHQuad_internal(p)}. Gradient (vector of
#'   derivatives with respect to each parameter) of \code{logLikNoAGHQuad(p)}.
#'   This is obtained using automatic differentiation (AD) with single-taping.
#'   First call will always be slower than later calls.
#'
#'   \item \code{gr_logLikNoAGHQuad(p)}. Gradient (vector of derivatives with
#'   respect to each parameter) of \code{logLikNoAGHQuad(p)}. This is obtained
#'   using automatic differentiation (AD) with double-taping. Results should
#'   match \code{gr_logLikNoAGHQuad_internal(p)} but may be more efficient after
#'   the first call.
#'
#'   \item \code{gr_AGHQuad(p)}. Gradient of the AGHQuad-approximated marginal
#'   likelihood at parameter value \code{p}.
#'
#'   \item \code{p_transformed_AGHQuad(pTransform)}. AGHQuad approximation at
#'         transformed (unconstrained) parameter value \code{pTransform}. To
#'         make maximizing the AGHQuad likelihood unconstrained, an automated
#'         transformation via \code{\link{parameterTransform}} is performed on
#'         any parameters with constraints indicated by their priors (even
#'         though the prior probabilities are not used).
#'
#' }
#'
#' Finally, methods that are primarily for internal use by other methods include:
#'
#' \itemize{
#'
#'    \item \code{p_transformed_gr_AGHQuad(pTransform)}. Gradient of the AGHQuad
#'     approximation (\code{p_transformed_AGHQuad(pTransform)}) at transformed 
#'     (unconstrained) parameter value \code{pTransform}.
#'
#'    \item \code{pInverseTransform(pTransform)}. Back-transform the transformed
#'    parameter value \code{pTransform} to original scale.
#'
#'    \item \code{derivspInverseTransform(pTransform, order)}. Derivatives of
#'    the back-transformation (i.e. inverse of parameter transformation) with
#'    respect to transformed parameters at \code{pTransform}. Derivative order
#'    is given by \code{order} (any of 0, 1, and/or 2).
#'
#'    \item \code{reInverseTransform(reTrans)}. Back-transform the transformed
#'    random effects value \code{reTrans} to original scale.
#'
#'    \item \code{derivsreInverseTransform(reTrans, order)}. Derivatives of the
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
#'    \item \code{inverseNegHess(p, reTransform)}. Calculate the inverse of the
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
#' }
#'
#' @section \code{control} list:
#' 
#' \code{buildAGHQuad} accepts the following control list elements:
#'
#' \itemize{
#'
#'   \item \code{split}. If TRUE (default), \code{randomEffectsNodes} will be
#'         split into conditionally independent sets if possible. This
#'         facilitates more efficient AGHQuad approximation because each
#'         conditionally independent set can be marginalized independently. If
#'         FALSE, \code{randomEffectsNodes} will be handled as one multivariate
#'         block, with one multivariate AGHQuad approximation. If \code{split}
#'         is a numeric vector, \code{randomEffectsNodes} will be split by
#'         \code{split}(\code{randomEffectsNodes}, \code{control$split}). The
#'         last option allows arbitrary control over how
#'         \code{randomEffectsNodes} are blocked.
#'
#'   \item \code{warn}. If TRUE (default), a warning is issued if
#'   \code{randomEffectsNodes} and/or \code{calcNodes} are provided but have
#'   extra or missing elements relative to their defaults.
#'
#'   \item \code{innerOptimControl}. See \code{optimControl}.
#'
#'   \item \code{innerOptimMethod}. See \code{optimMethod}.
#'
#'   \item \code{innerOptimStart}. see \code{optimStart}.
#'
#'   \item \code{outOptimControl}. A list of control parameters for maximizing
#'         the AGHQuad log-likelihood using \code{optim}. See 'Details' of
#'         \code{\link{optim}} for further information.
#'
#'   \item \code{allowNonPriors} If FALSE (default), all parameter nodes must
#'   have stochastic declarations (i.e. priors). These are not calculated as
#'   part of AGHQuad approximation, but they are used to determine boundaries on
#'   the valid ranges of parameters. If TRUE, parameters do not need
#'   declarations, and those without one will be assumed to have no constraints
#'   on valid values.
#'
#' }
#'
#' @section \code{AGHQuad_BASE}:
#' 
#' AGHQuad base class, upon which specific AGHQuad algorithm classes are based
#' by including \code{contains = AGHQuad_BASE}. This declares a list of
#' nimbleFunctions for a single AGHQuad approximation. Intended for internal use
#' only.
#'
#' @section \code{nimOneAGHQuad1D}:
#' 
#' This function constructs a single AGHQuad approximation when
#' \code{randomEffectsNodes} contains only one scalar node. It is mostly for
#' internal use by `buildAGHQuad`, when a scalar random effect is conditionally
#' independent from any other random effects. To use it directly, one has to
#' provide full inputs for all the arguments; there are no defaults.
#' 
#' This function generates an object that comprises a set of methods (functions), 
#' each accomplishing one piece of many calculations to obtain the AGHQuad 
#' approximation and its gradient w.r.t. model parameters. 
#'
#' Among these methods, six are most useful to a user:
#'
#' \itemize{
#'
#'   \item \code{AGHQuad1(p)}. AGHQuad approximation evaluated at the parameter
#'   value \code{p}. This function uses single AD taping for gradient and
#'   Hessian calculations and separate components.
#'
#'   \item \code{AGHQuad2(p)}. AGHQuad approximation evaluated at the parameter
#'   value \code{p}. This function uses double AD taping for gradient and
#'   Hessian calculations and separate components.
#'
#'   \item \code{AGHQuad3(p)}. AGHQuad approximation evaluated at the parameter
#'   value \code{p}. This function uses double AD taping for gradient and
#'   Hessian calculations and packs everything together.
#'
#'   \item \code{gr_AGHQuad1(p)}. Gradient of \code{AGHQuad1} with respect to
#'   parameters evaluated at the parameter value \code{p}.
#'
#'   \item \code{gr_AGHQuad2(p)}. Gradient of \code{AGHQuad2} with respect to
#'   parameters evaluated at the parameter value \code{p}.
#'
#'   \item \code{gr_AGHQuad3(p)}. Gradient of \code{AGHQuad3} with respect to
#'   parameters evaluated at the parameter value \code{p}.
#'
#' }
#' 
#' @section \code{nimOneAGHQuad}:
#' 
#' This function constructs a single AGHQuad approximation when
#' \code{randomEffectsNodes} contains more than one dimension (i.e. has a
#' non-scalar node and/or multiple scalar nodes). It is mostly for internal use
#' by `buildAGHQuad`. To use it directly, one has to provide full inputs for all
#' the arguments; there are no defaults.
#' 
#' The methods generated by this function are the same as \code{nimOneAGHQuad1D}. 
#'
#' @author Paul van Dam-Bates, Perry de Valpine
#' 
#' @name AGHQuad
#' 
#' @aliases AGHQuad buildAGHQuad
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
#' # Build AGHQuad approximation
#' pumpAGHQuad <- buildAGHQuad(pump)
#' 
#' \dontrun{
#' # Compile the model
#' Cpump <- compileNimble(pump)
#' CpumpAGHQuad <- compileNimble(pumpAGHQuad, project = pump)
#' # Calculate MLEs of parameters
#' MLEres <- CpumpAGHQuad$AGHQuadMLE()
#' # Calculate estimates and standard errors for parameters and random effects on original scale
#' allres <- CpumpAGHQuad$summary(MLEres, calcRandomEffectsStdError = TRUE)
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
