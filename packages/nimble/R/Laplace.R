#' NIMBLE Laplace approximation
#'
#' Laplace base class
Laplace_BASE <- nimbleFunctionVirtual(
  run = function() {},
  methods = list(
    Laplace1 = function(p = double(1)){
      returnType(double())
    },
    gr_Laplace1 = function(p = double(1)){
      returnType(double(1))
    },
    Laplace2 = function(p = double(1)){
      returnType(double())
    },
    gr_Laplace2 = function(p = double(1)){
      returnType(double(1))
    },
    Laplace3 = function(p = double(1)){
      returnType(double())
    },
    gr_Laplace3 = function(p = double(1)){
      returnType(double(1))
    }
  )
)

# To-do
# 1. Cache max inner logLik result - Done
# 1b. Use previous opt params to init max_inner_logLik - Done
# 2. Use conditionallyIndependentSets - Done
# 3. Insert checks for one_time_fixes - Done
# 4. Make nimOneLaplace1D case - Done
# 5. Get nre from paramTrans - Done
# 6. Use reInit name - Done
# 7. Group terms in %*% %*% - Done
# 8. Set up reset option for propagation, including for caches
# 9. Draft roxygen entry
# 10. p_transformed_gr needs gradient of inverse transformation - Done
# 11. Make all atomics use new - Done
# 12. test on spatial model
# 13. Make unified interface

nimOneLaplace1D <- nimbleFunction(
  contains = Laplace_BASE,
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes) {
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
    
    ## Automated transformation for random effects to ensure range of (-Inf, Inf) for each.
    reTrans <- parameterTransform(model, randomEffectsNodes)
    
    ## Numbers of random effects and parameters
    npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))
    nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))

    if(length(nre) != 1)
      stop("Number of random effects for nimOneLaplace1D must be 1.")
    
    ## Derivatives w.r.t. paramNodes and randomEffectsNodes are needed
    wrtNodes <- c(paramNodes, randomEffectsNodes)

    ## Indices of randomEffectsNodes inside wrtNodes
    re_indices <- as.numeric(c(npar+1, -1))
    
    ## Indices of paramNodes inside wrtNodes
    if(npar > 1) p_indices <- as.numeric(1:npar)
    else p_indices <- as.numeric(c(1, -1))
    
    ## Indices of randomEffectsNodes inside randomEffectsNodes for use in getting the derivative of 
    ## the inner log-likelihood (paramNodes fixed) w.r.t. randomEffectsNodes. 
    re_indices_inner <- as.numeric(c(1, -1))

    p_and_re_indices <- as.numeric(1:(npar + 1))
    ## Optim control list 
    optimControlList <- nimOptimDefaultControl()
  
    inner_makeUpdateNodes <- nimble:::makeUpdateNodes(wrtNodes = randomEffectsNodes, calcNodes = calcNodes, model = model)
    inner_updateNodes     <- inner_makeUpdateNodes$updateNodes
    inner_constantNodes   <- inner_makeUpdateNodes$constantNodes

    joint_makeUpdateNodes <- nimble:::makeUpdateNodes(wrtNodes = wrtNodes, calcNodes = calcNodes, model = model)
    joint_updateNodes     <- joint_makeUpdateNodes$updateNodes
    joint_constantNodes   <- joint_makeUpdateNodes$constantNodes

    ## The following are a set of member data
    ## for Laplace and gr_Laplace to store intermediate values
    ## for purposes of checking and comparing in case of bugs
    ## These are only used if record_intermediates_for_checking is set to TRUE
    record_intermediates_for_checking <- FALSE
    Lap_logdetNegHess_s <- numeric(1)
    Lap_opt_val_s <- numeric(1)
    gr_Lap_negHessian_s <- matrix(rnorm(4), nrow = 2)
    gr_Lap_grlogdetNegHesswrtp_s <- numeric(2)
    gr_Lap_grlogdetNegHesswrtre_s <- numeric(2)
    gr_Lap_hesslogLikwrtpre_s <- matrix(rnorm(4), nrow = 2)
    gr_Lap_gr_joint_logLik_wrt_p_s <- numeric(2)

    ## The following are used for caching values and gradient in the Laplace3 system
    Laplace3_saved_value <- numeric(1)
    Laplace3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    Laplace3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))

    ## The following are used for caching values and gradient in the Laplace3 system
    max_inner_logLik_saved_par <- as.numeric(c(1, -1))
    max_inner_logLik_saved_value <- numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    
    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE
    innerMethod <- "CG"
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
      return(max_inner_logLik_saved_par)
      returnType(double(1))
    },
    one_time_fixes = function() {
      ## Run this once after compiling; remove extraneous -1 if necessary
      if(one_time_fixes_done) return()
      re_indices <<- fix_one_vec(re_indices)
      re_indices_inner <<- fix_one_vec(re_indices_inner)
      max_inner_logLik_saved_par <<- fix_one_vec(max_inner_logLik_saved_par)
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        Laplace3_saved_gr <<- fix_one_vec(Laplace3_saved_gr)
        Laplace3_previous_p <<- fix_one_vec(Laplace3_previous_p)
        max_inner_logLik_previous_p <<- fix_one_vec(max_inner_logLik_previous_p)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit( reInit )
      one_time_fixes_done <<- TRUE
    },
    ## Joint log-likelihood with values of parameters fixed: used only for inner optimization
    inner_logLik = function(reTransform = double(1)) {
      re <- reTrans$inverseTransform(reTransform)
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(calcNodes) + reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    # Gradient of the joint log-likelihood (p fixed) w.r.t. transformed random effects: used only for inner optimization
    gr_inner_logLik_internal = function(reTransform = double(1)) {
      ans <- derivs(inner_logLik(reTransform), wrt = re_indices_inner, order = 1, model = model,
                    updateNodes   = inner_updateNodes, 
                    constantNodes = inner_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    ## WZ: it seems that the double-taping derivative methods sometimes give incorrect results.
    ## Currently derivative methods without double taping are used below.
    gr_inner_logLik = function(reTransform = double(1)) {
      ans <- derivs(gr_inner_logLik_internal(reTransform), wrt = re_indices_inner, order = 0, model = model,
                    updateNodes   = inner_updateNodes,
                    constantNodes = inner_updateNodes)
      return(ans$value)
      returnType(double(1))
    },
    # Setting optim control list as below does not work
    # ## Set up optim control list
    # ## Run this method to set up a new optim control list
    # set_optim_control = function(newControl = optimControlNimbleList()){
    #   optimControlList <<- newControl
    # },
    # ## Get optim control list
    # get_optim_control = function(){
    #   return(optimControlList)
    #   returnType(optimControlNimbleList())
    # },
    ## Solve the inner optimization for Laplace approximation
    max_inner_logLik = function(p = double(1)) {
      values(model, paramNodes) <<- p
      ## Random effects values from the model
      reInitTrans <- get_reInitTrans()
      ## control <- get_optim_control()
      optimControl <- nimOptimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 100 # 5000
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik, method = "CG", control = optimControl)
      ## Need to consider the case where optim does not converge
      if(optRes$convergence != 0) {
        print("Warning: optim does not converge for the inner optimization of Laplace approximation")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    },
    max_inner_logLik_internal = function(p = double(1)) {
      values(model, paramNodes) <<- p
      ## Random effects values from the model
      ## re <- values(model, randomEffectsNodes)
      ## reTransform <- reTrans$transform(re)
      reInitTrans <- get_reInitTrans()
      ## control <- get_optim_control()
      optimControl <- nimOptimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 5000
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik_internal, method = innerMethod, control = optimControl)
      ## Need to consider the case where optim does not converge
      if(optRes$convergence != 0) {
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
    },
    update_max_inner_logLik_internal = function(p = double(1)) {
      optRes <- max_inner_logLik_internal(p)
      max_inner_logLik_saved_par <<- optRes$par
      max_inner_logLik_saved_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
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
                    updateNodes   = joint_updateNodes, 
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = p_indices, order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_updateNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 1st order partial derivative w.r.t. transformed random effects
    gr_joint_logLik_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes   = joint_updateNodes, 
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 2nd order mixed partial derivative w.r.t. parameters and transformed random effects
    hess_joint_logLik_wrt_p_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    hess_joint_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      derivmat <- matrix(value = ans$value, nrow = npar) 
      return(derivmat)
      returnType(double(2))
    },
    ## Negative Hessian: 2nd order unmixed partial derivative w.r.t. transformed random effects
    negHess_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      ##hessian <- ans$jacobian
      return(-ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    negHess = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(negHess_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      neghess <- matrix(ans$value, nrow = nre)
      return(neghess)
      returnType(double(2))
    },
    ## Logdet negative Hessian
    logdetNegHess = function(p = double(1), reTransform = double(1)) {
      ## res <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = re_indices, order = 1, model = model,
      ##               updateNodes   = joint_updateNodes, 
      ##               constantNodes = joint_constantNodes)
      ## negHessian <- -(res$jacobian)
      negHessian <- negHess(p, reTransform)
      ans <- log(negHessian[1,1])
      ## ans <- logdet(negHess) ## Equivalent
      return(ans)
      returnType(double())
    },
    ## Gradient of logdet (negative) Hessian w.r.t. parameters
    gr_logdetNegHess_wrt_p_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = p_indices, order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_p_internal(p, reTransform), wrt = p_indices, order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Gradient of logdet (negative) Hessian w.r.t. transformed random effects
    gr_logdetNegHess_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = re_indices, order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_re_internal(p, reTransform), wrt = re_indices, order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ##
    joint_logLik_with_grad_and_hess_test = function(p = double(1), reTransform = double(1)) {
      # This returns a vector of  concatenated key quantities (see comment below for details) 
      # reTransform is the arg max of the inner logLik
      # We could consider returning only upper triangular elements of chol(-Hessian),
      #  and re-constituting as a matrix when needed.
      joint_logLik_res <- derivs(joint_logLik(p, reTransform),
                                 wrt = p_and_re_indices, # 1:(npar + nre)
                                 order = c(1, 2),
                                 model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      
      returnType(ADNimbleList())
      return(joint_logLik_res)
    },
    joint_logLik_with_grad_and_hess = function(p = double(1), reTransform = double(1)) {
      # This returns a vector of  concatenated key quantities (see comment below for details) 
      # reTransform is the arg max of the inner logLik
      # We could consider returning only upper triangular elements of chol(-Hessian),
      #  and re-constituting as a matrix when needed.
      joint_logLik_res <- derivs(joint_logLik(p, reTransform),
                                 wrt = p_and_re_indices, # 1:(npar + nre)
                                 order = c(1, 2),
                                 model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      negHessValue <- -joint_logLik_res$hessian[npar + 1, npar + 1, 1]
      logdetNegHessAns <- log(negHessValue)
      hess_wrt_p_wrt_re <- matrix(init = FALSE, nrow = npar, ncol = nre)
      for(i in 1:npar)
##        for(j in 1:nre)
        hess_wrt_p_wrt_re[i, 1] <- joint_logLik_res$hessian[i, npar + 1, 1]
      

      
      ans <- c(joint_logLik_res$jacobian[1, 1:npar], # I experimented with swapped order.  It's not being first.  It's the actual alpha parameter that has a problem, and which just happens to be first
               logdetNegHessAns,
               negHessValue, # cholNegHess
               hess_wrt_p_wrt_re)
      ## Indices to components of this are:
      ## gr_joint_logLik_wrt_p = (1:npar)                    [size = npar]
      ## logdetNegHess         = npar + 1                    [1]
      ## cholNegHess           = npar + 1 + (1 : nre*nre)    [nre x nre]
      ## hess_wrt_p_wrt_re     = npar + 1 + nre*nre + (1:npar*nre)  [npar x nre]
      return(ans)
      returnType(double(1))
      # return a concatenated vector
    },
    joint_logLik_with_higher_derivs_test = function(p = double(1), reTransform = double(1)) {
      higher_order_deriv_res <- derivs(joint_logLik_with_grad_and_hess(p, reTransform),
                                       wrt = p_and_re_indices,
                                       order = c(0, 1),
                                       model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      # value gives results from joint_logLik_with_grad_and_hess
      # jacobian gives derivs of these outputs wrt (p, re).
      # We only need gradient of logdetNegHess, which is the
      #   (1 + npar + 1, given in that order for sanity) row of jacobian
      # Other rows of the jacobian are wasted, but when this function
      # is meta-taped and optimized (part of CppAD), those calculations should be omitted
      returnType(ADNimbleList())
      return(higher_order_deriv_res)
    },
    joint_logLik_with_higher_derivs = function(p = double(1), reTransform = double(1)) {
      higher_order_deriv_res <- derivs(joint_logLik_with_grad_and_hess(p, reTransform),
                                       wrt = p_and_re_indices,
                                       order = c(0, 1),
                                       model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
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
    update_Laplace3_with_gr_test = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik_with_higher_derivs(p, reTransform), wrt = p_and_re_indices, order = 0,
                    model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans)
      returnType(ADNimbleList())
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
      #chol_negHess <- matrix(ans$value[(ind):(ind + nre*nre - 1)], nrow = nre, ncol = nre)
      negHessValue <- ans$value[ind]
      ind <- ind + 1
      hess_cross_terms <- numeric(value = ans$value[(ind):(ind + npar*1 - 1)], length = npar*1)
      ind <- ind + npar*1
      gr_logdetNegHess_wrt_p_v <- numeric(value = ans$value[(ind):(ind + npar - 1)], length = npar)
      ind <- ind + npar
      gr_logdetNegHess_wrt_re_v <- ans$value[ind]

      Laplace_value <- maxValue - 0.5 * logdetNegHess_value + 0.5 * 1 * log(2*pi)
      Laplace3_saved_value <<- Laplace_value 
      # print(Laplace_value)
      
      gr_Laplace_v <- gr_logLik_wrt_p - 0.5*(gr_logdetNegHess_wrt_p_v +  hess_cross_terms * (gr_logdetNegHess_wrt_re_v / negHessValue) )
      ## print( gr_logLik_wrt_p )
      ## print( gr_logdetNegHess_wrt_p_v )
      ## print( hess_cross_terms)
      ## print( gr_logdetNegHess_wrt_re_v)
      ## print( negHessValue)
      ## print( gr_Laplace_v )
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
      return(Laplace3_saved_value)
      returnType(double())
    },
    gr_Laplace3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      Laplace3_update(p)
      return(Laplace3_saved_gr)
      returnType(double(1))
    },
    ## Laplace approximation
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
      ## if(record_intermediates_for_checking) {
      ##   Lap_logdetNegHess_s <<- logdetNegHessian
      ##   Lap_opt_val_s <<- maxValue
      ## }
      return(ans)
      returnType(double())
    },
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
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation w.r.t. parameters
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
      ## if(record_intermediates_for_checking) {
      ##   gr_Lap_negHessian_s <<- negHessian
      ##   gr_Lap_grlogdetNegHesswrtp_s <<- grlogdetNegHesswrtp
      ##   gr_Lap_grlogdetNegHesswrtre_s <<- grlogdetNegHesswrtre
      ##   gr_Lap_hesslogLikwrtpre_s <<- hesslogLikwrtpre
      ##   gr_Lap_gr_joint_logLik_wrt_p_s <<- gr_joint_logLik_wrt_p(p, reTransform)
      ## }
      p1 <- gr_joint_logLik_wrt_p(p, reTransform) 
      ans <- p1 - 
        0.5 * (grlogdetNegHesswrtp + hesslogLikwrtpre * (grlogdetNegHesswrtre / negHessian))
      ## print( p1 )
      ## print( grlogdetNegHesswrtp )
      ## print( hesslogLikwrtpre)
      ## print( grlogdetNegHesswrtre)
      ## print( negHessian)
      ## print(ans )
      return(ans)
      returnType(double(1))
    },
    ## Gradient of the Laplace approximation w.r.t. parameters
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
  enableDerivs = list(inner_logLik                      = list(),
                      joint_logLik                      = list(),
                      gr_joint_logLik_wrt_re            = list(),
                      negHess                           = list(),
                      logdetNegHess                     = list(noDeriv_vars = c("d","i")),
                      gr_inner_logLik_internal          = list(),
                      gr_joint_logLik_wrt_p_internal    = list(),
                      gr_joint_logLik_wrt_re_internal   = list(),
                      hess_joint_logLik_wrt_p_wrt_re_internal = list(),
                      negHess_internal                  = list(),
                      gr_logdetNegHess_wrt_p_internal   = list(),
                      gr_logdetNegHess_wrt_re_internal  = list(),
                      joint_logLik_with_grad_and_hess = list(noDeriv_vars = c("i","j")),
                      joint_logLik_with_higher_derivs = list())
) ## End of nimOneLaplace1D


#' A nimbleFunction that performs a single Laplace approximation
nimOneLaplace <- nimbleFunction(
  contains = Laplace_BASE,
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes) {
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
    
    ## Automated transformation for random effects to ensure range of (-Inf, Inf) for each.
    reTrans <- parameterTransform(model, randomEffectsNodes)
    
    ## Numbers of random effects and parameters
    npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))
    nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))

    nreTrans <- reTrans$getTransformedLength()
    ## Derivatives w.r.t. paramNodes and randomEffectsNodes are needed
    wrtNodes <- c(paramNodes, randomEffectsNodes)

    ## Indices of randomEffectsNodes inside wrtNodes
    if(nreTrans > 1) reTrans_indices <- as.numeric((npar+1):(npar+nreTrans))
    else reTrans_indices <- as.numeric(c(npar+1, -1)) ## Ensure this is a vector
    
    ## Indices of paramNodes inside wrtNodes
    if(npar > 1) p_indices <- as.numeric(1:npar)
    else p_indices <- as.numeric(c(1, -1))
    
    ## Indices of randomEffectsNodes inside randomEffectsNodes for use in getting the derivative of 
    ## the inner log-likelihood (paramNodes fixed) w.r.t. randomEffectsNodes. 
    if(nreTrans > 1) reTrans_indices_inner <- as.numeric(1:nreTrans) 
    else reTrans_indices_inner <- as.numeric(c(1, -1))

    p_and_reTrans_indices <- as.numeric(1:(npar + nreTrans))
    ## Optim control list 
    optimControlList <- nimOptimDefaultControl()
  
    inner_makeUpdateNodes <- nimble:::makeUpdateNodes(wrtNodes = randomEffectsNodes, calcNodes = calcNodes, model = model)
    inner_updateNodes     <- inner_makeUpdateNodes$updateNodes
    inner_constantNodes   <- inner_makeUpdateNodes$constantNodes

    joint_makeUpdateNodes <- nimble:::makeUpdateNodes(wrtNodes = wrtNodes, calcNodes = calcNodes, model = model)
    joint_updateNodes     <- joint_makeUpdateNodes$updateNodes
    joint_constantNodes   <- joint_makeUpdateNodes$constantNodes

    ## The following are a set of member data
    ## for Laplace and gr_Laplace to store intermediate values
    ## for purposes of checking and comparing in case of bugs
    ## These are only used if record_intermediates_for_checking is set to TRUE
    ## record_intermediates_for_checking <- FALSE
    ## Lap_logdetNegHess_s <- numeric(1)
    ## Lap_opt_val_s <- numeric(1)
    ## gr_Lap_negHessian_s <- matrix(rnorm(4), nrow = 2)
    ## gr_Lap_grlogdetNegHesswrtp_s <- numeric(2)
    ## gr_Lap_grlogdetNegHesswrtre_s <- numeric(2)
    ## gr_Lap_hesslogLikwrtpre_s <- matrix(rnorm(4), nrow = 2)
    ## gr_Lap_gr_joint_logLik_wrt_p_s <- numeric(2)

    ## The following are used for caching values and gradient in the Laplace3 system
    Laplace3_saved_value <- numeric(1)
    Laplace3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    Laplace3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))

    ## The following are used for caching values and gradient in the Laplace3 system
    max_inner_logLik_saved_par <- if(nreTrans > 1) numeric(nreTrans) else as.numeric(c(1, -1))
    max_inner_logLik_saved_value <- numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    
    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE
    num_times <- 11
    times <- rep(0, num_times)
    num_inner_times <- 1000
    i_inner_time <- 1
    inner_logLik_times <- rep(0, num_inner_times)
    innerMethod <- "BFGS"
  },
  run = function(){},
  methods = list(
    get_times = function() {return(times); returnType(double(1))},
    reset_times = function() {times <<- rep(0, num_times)},
    add_time = function(t = double(), i = integer()) {times[i] <<- times[i] + t},
    
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
      return(max_inner_logLik_saved_par)
      returnType(double(1))
    },
    one_time_fixes = function() {
      ## Run this once after compiling; remove extraneous -1 if necessary
      if(one_time_fixes_done) return()
      if(nre == 1) {
        reTrans_indices <<- fix_one_vec(reTrans_indices)
        reTrans_indices_inner <<- fix_one_vec(reTrans_indices_inner)
        max_inner_logLik_saved_par <<- fix_one_vec(max_inner_logLik_saved_par)
      }
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        Laplace3_saved_gr <<- fix_one_vec(Laplace3_saved_gr)
        Laplace3_previous_p <<- fix_one_vec(Laplace3_previous_p)
        max_inner_logLik_previous_p <<- fix_one_vec(max_inner_logLik_previous_p)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit( reInit )
      one_time_fixes_done <<- TRUE
    },
    ## Joint log-likelihood with values of parameters fixed: used only for inner optimization
    inner_logLik = function(reTransform = double(1)) {
      re <- reTrans$inverseTransform(reTransform)
      values(model, randomEffectsNodes) <<- re
      ans <- model$calculate(calcNodes) + reTrans$logDetJacobian(reTransform)
      return(ans)
      returnType(double())
    },
    # Gradient of the joint log-likelihood (p fixed) w.r.t. transformed random effects: used only for inner optimization
    gr_inner_logLik_internal = function(reTransform = double(1)) {
      ans <- derivs(inner_logLik(reTransform), wrt = reTrans_indices_inner,
                    order = 1, model = model,
                    updateNodes   = inner_updateNodes, 
                    constantNodes = inner_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping for efficiency
    ## WZ: it seems that the double-taping derivative methods sometimes give incorrect results.
    ## Currently derivative methods without double taping are used below.
    gr_inner_logLik = function(reTransform = double(1)) {
      t10 <- run.time( {
      ans <- derivs(gr_inner_logLik_internal(reTransform), wrt = reTrans_indices_inner,
                    order = 0, model = model,
                    updateNodes   = inner_updateNodes,
                    constantNodes = inner_updateNodes)
      })
      add_time(t10, 10)
      return(ans$value)
      returnType(double(1))
    },
    # Setting optim control list as below does not work
    # ## Set up optim control list
    # ## Run this method to set up a new optim control list
    # set_optim_control = function(newControl = optimControlNimbleList()){
    #   optimControlList <<- newControl
    # },
    # ## Get optim control list
    # get_optim_control = function(){
    #   return(optimControlList)
    #   returnType(optimControlNimbleList())
    # },
    ## Solve the inner optimization for Laplace approximation
    max_inner_logLik = function(p = double(1)) {
      values(model, paramNodes) <<- p
      ## Random effects values from the model
      reInitTrans <- get_reInitTrans()
      
      # If optim can't get a value from initial parameters (random effects),
      # it will error out and we don't have a good way to catch it.
      # This might be possible via but I am not sure how.  It is thrown from a direct error() call in the optim C++ functions.
      # Hence we check the initial parameters (random effects) and return -Inf manually.
      # Note that what could be happening is the (actual) parameters from the outer optimization
      # could be invalid, which would mean that no values of initial parameters (random effects)
      # here in the inner optimization will be valid.  The outer optim can handle an Inf and search other parameter values,
      # but we need the inner optimization to exit gracefully if it can't do anything.
      fn_init <- inner_logLik(reInitTrans)
##      print("in max_inner_logLik\n")
      if((fn_init == Inf) | (fn_init == -Inf) | (is.nan(fn_init)) | (is.na(fn_init))) {
##        print("init params failed\n")
        optRes <- optimResultNimbleList$new()
        optRes$par <- reInitTrans
        optRes$value <- -Inf
        optRes$convergence <- -1
        return(optRes)
      }

      ## control <- get_optim_control()
      optimControl <- nimOptimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 5000
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik,
                      method = innerMethod, control = optimControl)
      ## Need to consider the case where optim does not converge
      if(optRes$convergence != 0) {
        print("Warning: optim does not converge for the inner optimization of Laplace approximation")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    },
    max_inner_logLik_internal = function(p = double(1)) {
      values(model, paramNodes) <<- p
      ## Random effects values from the model
      reInitTrans <- get_reInitTrans()

      ## control <- get_optim_control()
      optimControl <- nimOptimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 5000
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik_internal,
                      method = "CG", control = optimControl)
      ## Need to consider the case where optim does not converge
      if(optRes$convergence != 0) {
        print("Warning: optim does not converge for the inner optimization of Laplace approximation")
      }
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## These two update methods for max_inner_logLik use the same member data caches
    update_max_inner_logLik = function(p = double(1)) {
      t9 <- run.time( {
        optRes <- max_inner_logLik(p)
      })
      add_time(t9, 9)

      inner_logLik_times[i_inner_time] <<- t9
      if(i_inner_time == length(inner_logLik_times)) {
        inner_logLik_times <<- c(inner_logLik_times, rep(0, 1000))
      }
      i_inner_time <<- i_inner_time + 1
      
      max_inner_logLik_saved_par <<- optRes$par
      max_inner_logLik_saved_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
    },
    update_max_inner_logLik_internal = function(p = double(1)) {
      optRes <- max_inner_logLik_internal(p)
      max_inner_logLik_saved_par <<- optRes$par
      max_inner_logLik_saved_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
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
      ans <- derivs(joint_logLik(p, reTransform), wrt = p_indices,
                    order = 1, model = model,
                    updateNodes   = joint_updateNodes, 
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = p_indices,
                    order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_updateNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 1st order partial derivative w.r.t. transformed random effects
    gr_joint_logLik_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik(p, reTransform), wrt = reTrans_indices,
                    order = 1, model = model,
                    updateNodes   = joint_updateNodes, 
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_joint_logLik_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = reTrans_indices,
                    order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## 2nd order mixed partial derivative w.r.t. parameters and transformed random effects
    hess_joint_logLik_wrt_p_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_p_internal(p, reTransform), wrt = reTrans_indices,
                    order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    hess_joint_logLik_wrt_p_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform), wrt = reTrans_indices,
                    order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      derivmat <- matrix(value = ans$value, nrow = npar) 
      return(derivmat)
      returnType(double(2))
    },
    ## Negative Hessian: 2nd order unmixed partial derivative w.r.t. transformed random effects
    negHess_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_joint_logLik_wrt_re_internal(p, reTransform), wrt = reTrans_indices,
                    order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      ##hessian <- ans$jacobian
      return(-ans$jacobian)
      returnType(double(2))
    },
    ## Double taping
    negHess = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(negHess_internal(p, reTransform), wrt = reTrans_indices,
                    order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      neghess <- matrix(ans$value, nrow = nreTrans)
      return(neghess)
      returnType(double(2))
    },
    ## This was an experiment in double taping the second order.  It was very slow. 
    ## negHess2_internal = function(p = double(1), reTransform = double(1)) {
    ##   ans <- derivs(joint_logLik(p, reTransform), wrt = reTrans_indices,
    ##                 order = 2, model = model,
    ##                 updateNodes   = joint_updateNodes,
    ##                 constantNodes = joint_constantNodes)
    ##   neghess <- matrix(nrow = nreTrans, ncol = nreTrans)
    ##   for(i in 1:nreTrans)
    ##     for(j in 1:nreTrans)
    ##       neghess[i, j] <- -ans$hessian[i, j, 1]
    ##   return(neghess)
    ##   returnType(double(2))
    ## },
    ## negHess2 = function(p = double(1), reTransform = double(1)) {
    ##   ans <- derivs(negHess2_internal(p, reTransform), wrt = reTrans_indices,
    ##                 order = 0, model = model,
    ##                 updateNodes   = joint_updateNodes,
    ##                 constantNodes = joint_constantNodes)
    ##   neghess <- matrix(ans$value, nrow = nreTrans)
    ##   return(neghess)
    ##   returnType(double(2))
    ## },
    ## Logdet negative Hessian
    logdetNegHess = function(p = double(1), reTransform = double(1)) {
      negHessian <- negHess(p, reTransform)
      cholNegHess <- chol(negHessian)
      ans <- 2 * sum(log(diag(cholNegHess)))
      ## ans <- logdet(negHess) ## Equivalent
      return(ans)
      returnType(double())
    },
    ## Gradient of logdet (negative) Hessian w.r.t. parameters
    gr_logdetNegHess_wrt_p_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = p_indices,
                    order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_p = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_p_internal(p, reTransform), wrt = p_indices,
                    order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ## Gradient of logdet (negative) Hessian w.r.t. transformed random effects
    gr_logdetNegHess_wrt_re_internal = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(logdetNegHess(p, reTransform), wrt = reTrans_indices,
                    order = 1, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$jacobian[1,])
      returnType(double(1))
    },
    ## Double taping
    gr_logdetNegHess_wrt_re = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(gr_logdetNegHess_wrt_re_internal(p, reTransform), wrt = reTrans_indices,
                    order = 0, model = model,
                    updateNodes   = joint_updateNodes,
                    constantNodes = joint_constantNodes)
      return(ans$value)
      returnType(double(1))
    },
    ##
    joint_logLik_with_grad_and_hess_test = function(p = double(1), reTransform = double(1)) {
      # This returns a vector of  concatenated key quantities (see comment below for details) 
      # reTransform is the arg max of the inner logLik
      # We could consider returning only upper triangular elements of chol(-Hessian),
      #  and re-constituting as a matrix when needed.
      joint_logLik_res <- derivs(joint_logLik(p, reTransform),
                                 wrt = p_and_reTrans_indices, # 1:(npar + nreTrans)
                                 order = c(1, 2),
                                 model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      
      returnType(ADNimbleList())
      return(joint_logLik_res)
    },
    joint_logLik_with_grad_and_hess = function(p = double(1), reTransform = double(1)) {
      # This returns a vector of  concatenated key quantities (see comment below for details) 
      # reTransform is the arg max of the inner logLik
      # We could consider returning only upper triangular elements of chol(-Hessian),
      #  and re-constituting as a matrix when needed.
      joint_logLik_res <- derivs(joint_logLik(p, reTransform),
                                 wrt = p_and_reTrans_indices, # 1:(npar + nreTrans)
                                 order = c(1, 2),
                                 model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      negHessUpper <- matrix(init = FALSE, nrow = nre, ncol = nreTrans)
      for(i in 1:nreTrans)
        for(j in i:nreTrans)
          negHessUpper[i,j] <- -joint_logLik_res$hessian[npar + i, npar + j, 1]
      cholNegHess <- chol(negHessUpper)
      logdetNegHessAns <- 2 * sum(log(diag(cholNegHess)))
      hess_wrt_p_wrt_re <- matrix(init = FALSE, nrow = npar, ncol = nre)
      for(i in 1:npar)
        for(j in 1:nreTrans)
          hess_wrt_p_wrt_re[i, j] <- joint_logLik_res$hessian[i, npar + j, 1]
      
      ans <- c(joint_logLik_res$jacobian[1, 1:npar], # I experimented with swapped order.  It's not being first.  It's the actual alpha parameter that has a problem, and which just happens to be first
               logdetNegHessAns,
               cholNegHess,
               hess_wrt_p_wrt_re)
      ## Indices to components of this are:
      ## gr_joint_logLik_wrt_p = (1:npar)                    [size = npar]
      ## logdetNegHess         = npar + 1                    [1]
      ## cholNegHess           = npar + 1 + (1 : nreTrans * nreTrans)    [nreTrans x nreTrans]
      ## hess_wrt_p_wrt_re     = npar + 1 + nre*nre + (1:npar*nreTrans)  [npar x nreTrans]
      return(ans)
      returnType(double(1))
      # return a concatenated vector
    },
    joint_logLik_with_higher_derivs_test = function(p = double(1), reTransform = double(1)) {
      higher_order_deriv_res <- derivs(joint_logLik_with_grad_and_hess(p, reTransform),
                                       wrt = p_and_reTrans_indices,
                                       order = c(0, 1),
                                       model = model, updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      # value gives results from joint_logLik_with_grad_and_hess
      # jacobian gives derivs of these outputs wrt (p, re).
      # We only need gradient of logdetNegHess, which is the
      #   (1 + npar + 1, given in that order for sanity) row of jacobian
      # Other rows of the jacobian are wasted, but when this function
      # is meta-taped and optimized (part of CppAD), those calculations should be omitted
      returnType(ADNimbleList())
      return(higher_order_deriv_res)
    },
    joint_logLik_with_higher_derivs = function(p = double(1), reTransform = double(1)) {
      higher_order_deriv_res <- derivs(joint_logLik_with_grad_and_hess(p, reTransform),
                                       wrt = p_and_reTrans_indices,
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
    update_Laplace3_with_gr_test = function(p = double(1), reTransform = double(1)) {
      ans <- derivs(joint_logLik_with_higher_derivs(p, reTransform), wrt = p_and_reTrans_indices,
                    order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      return(ans)
      returnType(ADNimbleList())
    },
    update_Laplace3_with_gr = function(p = double(1), reset = logical(0, default = FALSE)) {
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value        

      ans <- derivs(joint_logLik_with_higher_derivs(p, reTransform), wrt = p_and_reTrans_indices,
                    order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      ind <- 1
      # all "logLik" here is joint log likelihood (i.e. for p and re)
      gr_logLik_wrt_p <- ans$value[(ind):(ind + npar - 1)]
      ind <- ind + npar
      logdetNegHess_value <- ans$value[ind]
      ind <- ind + 1
      chol_negHess <- matrix(ans$value[(ind):(ind + nreTrans*nreTrans - 1)],
                             nrow = nreTrans, ncol = nreTrans)
      ind <- ind + nreTrans*nreTrans
      hess_cross_terms <- matrix(ans$value[(ind):(ind + npar*nreTrans - 1)],
                                 nrow = npar, ncol = nreTrans)
      ind <- ind + npar*nreTrans
      gr_logdetNegHess_wrt_p_v <- ans$value[(ind):(ind + npar - 1)]
      ind <- ind + npar
      gr_logdetNegHess_wrt_re_v <- ans$value[(ind):(ind + nreTrans - 1)]

      Laplace_value <- maxValue - 0.5 * logdetNegHess_value + 0.5 * nreTrans * log(2*pi)
      Laplace3_saved_value <<- Laplace_value 
      # print(Laplace_value)
      
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
      return(Laplace3_saved_value)
      returnType(double())
    },
    gr_Laplace3 = function(p = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      Laplace3_update(p)
      return(Laplace3_saved_gr)
      returnType(double(1))
    },
    ## Laplace approximation
    Laplace2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value        

      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed

      # time 1
      t1 <- run.time( {
        logdetNegHessian <- logdetNegHess(p, reTransform)
      })
      add_time(t1, 1)
      ## Laplace approximation
      t2 <- run.time( {
        ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      })
      add_time(t2, 2)
      ## if(record_intermediates_for_checking) {
      ##   Lap_logdetNegHess_s <<- logdetNegHessian
      ##   Lap_opt_val_s <<- maxValue
      ## }
      return(ans)
      returnType(double())
    },
    Laplace1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_saved_par
      maxValue <- max_inner_logLik_saved_value        

      logdetNegHessian <- logdetNegHess(p, reTransform)
      ## Laplace approximation
      ans <- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation w.r.t. parameters
    gr_Laplace2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_saved_par

      t3 <- run.time( {
        negHessian <- negHess(p, reTransform)
      })
      add_time(t3, 3)
      t4 <- run.time( {
        invNegHessian <- inverse(negHessian)
      })
      add_time(t4, 4)
      t5 <- run.time( {
        grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p(p, reTransform)
      })
      add_time(t5, 5)
      t6 <- run.time( {
        grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re(p, reTransform)
      })
      add_time(t6, 6)
      t7 <- run.time( {
        hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re(p, reTransform)
      })
      add_time(t7, 7)
      ## if(record_intermediates_for_checking) {
      ##   gr_Lap_negHessian_s <<- negHessian
      ##   gr_Lap_grlogdetNegHesswrtp_s <<- grlogdetNegHesswrtp
      ##   gr_Lap_grlogdetNegHesswrtre_s <<- grlogdetNegHesswrtre
      ##   gr_Lap_hesslogLikwrtpre_s <<- hesslogLikwrtpre
      ##   gr_Lap_gr_joint_logLik_wrt_p_s <<- gr_joint_logLik_wrt_p(p, reTransform)
      ## }
      t8 <- run.time({
      ans <- gr_joint_logLik_wrt_p(p, reTransform) - 
        0.5 * (grlogdetNegHesswrtp + (grlogdetNegHesswrtre %*% invNegHessian) %*% t(hesslogLikwrtpre))
      })
      add_time(t8, 8)
      
      return(ans[1,])
      returnType(double(1))
    },
    ## Gradient of the Laplace approximation w.r.t. parameters
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
  enableDerivs = list(inner_logLik                      = list(),
                      joint_logLik                      = list(),
                      gr_joint_logLik_wrt_re            = list(),
                      negHess                           = list(),
                      logdetNegHess                     = list(noDeriv_vars = c("d","i")),
                      gr_inner_logLik_internal          = list(),
                      gr_joint_logLik_wrt_p_internal    = list(),
                      gr_joint_logLik_wrt_re_internal   = list(),
                      hess_joint_logLik_wrt_p_wrt_re_internal = list(),
                      negHess_internal                  = list(),
                      ## negHess2_internal                  = list(noDeriv_vars = c("i","j")),
                      gr_logdetNegHess_wrt_p_internal   = list(),
                      gr_logdetNegHess_wrt_re_internal  = list(),
                      joint_logLik_with_grad_and_hess = list(noDeriv_vars = c("i","j")),
                      joint_logLik_with_higher_derivs = list())
) ## End of nimOneLaplace


#' From a set of nodes, find out nodes that are dependent of tarNodes
getDependentNodes <- function(model, nodes, tarNodes){
  nodes <- model$expandNodeNames(nodes)
  tarNodes <- model$expandNodeNames(tarNodes)
  ## Nodes that depend on tarNodes 
  depNodes <- intersect(model$getDependencies(tarNodes), nodes)
  ## Nodes that tarNodes depend on
  for(node in setdiff(nodes, depNodes)){
    if(any(tarNodes %in% model$getDependencies(node))){
      depNodes <- c(depNodes, node)
    }
  }
  return(depNodes)
}

#' Find a maximum set of conditionally dependent nodes from a set of nodes
findOneMaxDependentNodeSet <- function(model, nodes){
  nodes <- model$expandNodeNames(nodes)
  ## Dependent nodes of the first node
  depNodes <- getDependentNodes(model, nodes, nodes[1])
  ## Dependent nodes of depNodes
  ddepNodes <- getDependentNodes(model, nodes, depNodes)
  ## Add new dependent nodes into depNodes if any
  while(any(setdiff(nodes, depNodes) %in% ddepNodes)){
    depNodes <- ddepNodes
    ddepNodes <- getDependentNodes(model, nodes, depNodes)
  }
  return(depNodes)
}

#' Divide a set of nodes into several conditionally independent groups
findAllDependentNodeSets <- function(model, nodes){
  nodes <- model$expandNodeNames(nodes)
  res <- list()
  res[[1]] <- findOneMaxDependentNodeSet(model, nodes)
  groupedNodes <- res[[1]] ## Grouped nodes 
  nNodes <- length(nodes)
  nLeftNodes <- nNodes - length(groupedNodes) ## Number of ungrouped nodes
  i <- 1
  while(nLeftNodes > 0){
    i <- i + 1
    res[[i]] <- findOneMaxDependentNodeSet(model, setdiff(nodes, groupedNodes))
    nLeftNodes <- nLeftNodes - length(res[[i]])
    groupedNodes <- c(groupedNodes, res[[i]])
  }
  return(res)
}

#' nimbleFunction for Laplace approximation
#'
#' There are currently three versions:
#' LaplaceMLE1 and associated "Laplace1" functions: single-tapes of component calculations
#' LaplaceMLE2 and associated "Laplace2" functions: double-tapes of component calculations
#' LaplaceMLE3 and associated "Laplace3" functions: double-tape of combined calculations
#' LaplaceMLE uses the one of these determined by methodID
#' methodID can be changed by set_method and queried by get_method
#' 
#' paramNodes: default to top nodes
#' randomEffectsNodes: default to latent nodes that depend on paramNodes
#'                     A warning is issued if randomEffectsNodes is provided
#'                     and has extra or missing elements
#' calcNodes: default to model$geteDependencies(randomEffectsNodes)
#'            A warning is issued if calcNodes is provided and
#'            has extra or missing elements
#'
#' There may be deterministic nodes between paramNodes and
#' randomEffectsNodes.  These will be included in calculations automatically.
#'
#' control list elements:
#' split: If TRUE, randomEffectsNodes will be split into conditionally independent sets.
#'        If FALSE, randomEffectsNodes will be handled in on multivariate block.
#'        If a vector, randomEffectsNodes will be split by split(randomEffectsNodes, control$split)
#'          The last option allows arbitrary control over how randomEffectsNodes are blocked.
#'
#' Examples
#' laplace <- buildLaplace(model) # should do everything automatically
#' Claplace <- compileNimble(laplace, project = model)
#' Claplace$LaplaceMLE()
buildLaplace <- nimbleFunction(
  setup = function(model, paramNodes, randomEffectsNodes, calcNodes,
                   control = list()) {
    if(is.null(control$split)) split <- TRUE else split <- control$split
    if(is.null(control$warn)) warn <- TRUE else warn <- control$warn
    laplace_nfl <- nimbleFunctionList(Laplace_BASE)

    paramProvided <- !missing(paramNodes)
    reProvided <- !missing(randomEffectsNodes)
    calcProvided <- !missing(calcNodes)

    if(!paramProvided) {
      paramNodes <- model$getNodeNames(topOnly = TRUE)
    } else {
      paramNodes <- model$expandNodeNames(paramNodes)
    }

    if((!reProvided) || warn) {
      paramDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
      reNodesDefault <- model$getNodeNames(latentOnly = TRUE)
      reNodesDefault <- intersect(reNodesDefault, paramDeps)
    }
    if(reProvided)
      randomEffectsNodes <- model$expandNodeNames(randomEffectsNodes)
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
      calcNodesDefault <- model$getDependencies(randomEffectsNodes)
    }
    if(calcProvided)
      calcNodes <- model$expandNodeNames(calcNodes)
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
      if(length(calcCheck)) {
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
    if(!calcProvided) {
      calcNodes <- calcNodesDefault
    }

    # If split == FALSE, do all random effects in one set
    # If split is integer (primarily for testing)
    if(isFALSE(split)) {
      laplace_nfl[[1]] <- nimOneLaplace(model, paramNodes, randomEffectsNodes, calcNodes)
    }
    else {
      ## 1. Split randomEffectsNodes into sets
      if(isTRUE(split))
        reSets <- model$getConditionallyIndependentSets(nodes = randomEffectsNodes,
                                                        givenNodes = setdiff(c(paramNodes, calcNodes),
                                                                             randomEffectsNodes))
      else if(is.numeric(split))
        reSets <- split(randomEffectsNodes, split)
      else stop("Invalid value for \'split\' provided in control list")
      
      num_reSets <- length(reSets)
      if(num_reSets == 0)
        stop("There was a problem determining conditionally independent sets for this model.")
      for(i in seq_along(reSets)) {
        ## Work with one conditionally independent set of latent states
        these_reNodes <- reSets[[i]]
        ## find paramNodes and calcNodes for this set of reNodes
        these_reDeps <- model$getDependencies(these_reNodes) ## candidate calcNodes via reNodes
        these_calcNodes <- intersect(calcNodes, these_reDeps) ## definite calcNodes

        ## paramNodes are the same for all laplace_nfl elements.
        ## In the future this could be customized.
        if(length(model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE)) > 1) {
          laplace_nfl[[i]] <- nimOneLaplace(model, paramNodes,
                                            these_reNodes, these_calcNodes)
        } else {
          laplace_nfl[[i]] <- nimOneLaplace1D(model, paramNodes,
                                              these_reNodes, these_calcNodes)
        }
      }
    }
      
    npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))
    ## Automated transformation for parameters
    paramsTransformation <- parameterTransform(model, paramNodes)
    pTransform_length <- paramsTransformation$getTransformedLength()
    pTransform_indices <- if(pTransform_length > 1) 1:pTransform_length else c(1, -1)
    one_time_fixes_done <- FALSE

    methodID <- 3
  },
  run = function(){},
  methods = list(
    set_method = function(method = integer()) {
      methodID <<- method
    },
    get_method = function() {
      return(methodID)
      returnType(integer())
    },
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(pTransform_length == 1) {
        if(length(pTransform_indices) == 2)
          pTransform_indices <<- numeric(length = 1, value = 1)
      }
      one_time_fixes_done <<- TRUE
    },
    Laplace = function(p = double(1)) {
      ans <- -Inf
      if(methodID == 1)
        ans <- Laplace1(p)
      else if(methodID == 2)
        ans <- Laplace2(p)
      else if(methodID == 3)
        ans <- Laplace3(p)
      return(ans)
      returnType(double())
    },
    ## Laplace approximation
    Laplace1 = function(p = double(1)){
      ans <- 0
      for(i in seq_along(laplace_nfl)){
        ans <- ans + laplace_nfl[[i]]$Laplace1(p)
      }
      return(ans)
      returnType(double())
    },
    Laplace2 = function(p = double(1)){
      ans <- 0
      for(i in seq_along(laplace_nfl)){
        ans <- ans + laplace_nfl[[i]]$Laplace2(p)
      }
      return(ans)
      returnType(double())
    },
    Laplace3 = function(p = double(1)){
      ans <- 0
      for(i in seq_along(laplace_nfl)){
        ans <- ans + laplace_nfl[[i]]$Laplace3(p)
      }
      return(ans)
      returnType(double())
    },
    ## Gradient of the Laplace approximation w.r.t. parameters
    gr_Laplace = function(p = double(1)) {
      if(methodID == 1)
        ans <- gr_Laplace1(p)
      else if(methodID == 2)
        ans <- gr_Laplace2(p)
      else if(methodID == 3)
        ans <- gr_Laplace3(p)
      return(ans)
      returnType(double(1))
    },
    gr_Laplace1 = function(p = double(1)){
      ans <- numeric(length = npar)
      for(i in seq_along(laplace_nfl)){
        ans <- ans + laplace_nfl[[i]]$gr_Laplace1(p)
      }
      return(ans)
      returnType(double(1))
    },
    gr_Laplace2 = function(p = double(1)){
      ans <- numeric(length = npar)
      for(i in seq_along(laplace_nfl)){
        ans <- ans + laplace_nfl[[i]]$gr_Laplace2(p)
      }
      return(ans)
      returnType(double(1))
    },
    gr_Laplace3 = function(p = double(1)){
      ans <- numeric(length = npar)
      for(i in seq_along(laplace_nfl)){
        ans <- ans + laplace_nfl[[i]]$gr_Laplace3(p)
      }
      return(ans)
      returnType(double(1))
    },
    ## Laplace approximation in terms of transformed parameters
    p_transformed_Laplace = function(pTransform = double(1)) {
      if(methodID == 1)
        ans <- p_transformed_Laplace1(pTransform)
      else if(methodID == 2)
        ans <- p_transformed_Laplace2(pTransform)
      else if(methodID == 3)
        ans <- p_transformed_Laplace3(pTransform)
      return(ans)
      returnType(double())
    },
    p_transformed_Laplace1 = function(pTransform = double(1)) {
      p <- paramsTransformation$inverseTransform(pTransform)
      ans <- Laplace1(p)
      return(ans)
      returnType(double())
    },
    p_transformed_Laplace2 = function(pTransform = double(1)) {
      p <- paramsTransformation$inverseTransform(pTransform)
      ans <- Laplace2(p)
      if(is.nan(ans)) ans <- -Inf
      return(ans)
      returnType(double())
    },
    p_transformed_Laplace3 = function(pTransform = double(1)) {
      p <- paramsTransformation$inverseTransform(pTransform)
      ans <- Laplace3(p)
      return(ans)
      returnType(double())
    },

    inverseTransform = function(pTransform = double(1)) {
      p <- paramsTransformation$inverseTransform(pTransform)
      return(p)
      returnType(double(1))
    },
    derivsInverseTransform = function(pTransform = double(1),
                                      order = double(1)) {
      if(!one_time_fixes_done) one_time_fixes()
      ans <- derivs(inverseTransform(pTransform),
                    wrt = pTransform_indices,
                    order = order)
      return(ans)
      returnType(ADNimbleList())
    },
    
    ## Gradient of the Laplace approximation in terms of transformed parameters
    p_transformed_gr_Laplace = function(pTransform = double(1)) {
      if(methodID == 1)
        ans <- p_transformed_gr_Laplace1(pTransform)
      else if(methodID == 2)
        ans <- p_transformed_gr_Laplace2(pTransform)
      else if(methodID == 3)
        ans <- p_transformed_gr_Laplace3(pTransform)
      return(ans)
      returnType(double(1))
    },
    p_transformed_gr_Laplace1 = function(pTransform = double(1)) {
#      p <- paramsTransformation$inverseTransform(pTransform)
      pDerivs <- derivsInverseTransform(pTransform, c(0, 1))
      ans <- gr_Laplace1(pDerivs$value)
      ans <- (ans %*% pDerivs$jacobian)[1,]
      return(ans)
      returnType(double(1))
    },
    p_transformed_gr_Laplace2 = function(pTransform = double(1)) {
      pDerivs <- derivsInverseTransform(pTransform, c(0, 1))
      ans <- gr_Laplace2(pDerivs$value)
      ans <- (ans %*% pDerivs$jacobian)[1,]
      return(ans)
      returnType(double(1))
    },
    p_transformed_gr_Laplace3 = function(pTransform = double(1)) {
      pDerivs <- derivsInverseTransform(pTransform, c(0, 1))
      ans <- gr_Laplace3(pDerivs$value)
      ans <- (ans %*% pDerivs$jacobian)[1,]
      return(ans)
      returnType(double(1))
    },

    ## Find Laplace MLEs of parameters
    ## Version 1: plain single-taping of component calculations
    LaplaceMLE = function(pStart = double(1)) {
      if(methodID == 1)
        ans <- LaplaceMLE1(pStart)
      else if(methodID == 2)
        ans <- LaplaceMLE2(pStart)
      else if(methodID == 3)
        ans <- LaplaceMLE3(pStart)
      return(ans)
      returnType(optimResultNimbleList())
    },
    LaplaceMLE1 = function(pStart = double(1)){
      pStartTransform <- paramsTransformation$transform(pStart)
      optimControl <- optimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 5000
      optRes <- optim(pStartTransform,
                      p_transformed_Laplace1, p_transformed_gr_Laplace1,
                      method = "BFGS", control = optimControl)
      ## Return MLEs on the original scale
      transformedMLEs <- optRes$par
      MLEs <- paramsTransformation$inverseTransform(transformedMLEs)
      optRes$par <- MLEs
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Find Laplace MLEs of parameters
    ## Version 2: Double-taping of component calculations
    LaplaceMLE2 = function(pStart = double(1)){
      pStartTransform <- paramsTransformation$transform(pStart)
      optimControl <- optimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 5000
      optRes <- optim(pStartTransform,
                      p_transformed_Laplace2, p_transformed_gr_Laplace2,
                      method = "BFGS", control = optimControl)
      ## Return MLEs on the original scale
      transformedMLEs <- optRes$par
      MLEs <- paramsTransformation$inverseTransform(transformedMLEs)
      optRes$par <- MLEs
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Find Laplace MLEs of parameters
    ## Version 3: Double-taping of calculations packed together to reduce redundancy.
    LaplaceMLE3 = function(pStart = double(1)){
      pStartTransform <- paramsTransformation$transform(pStart)
      optimControl <- optimDefaultControl()
      optimControl$fnscale <- -1
      optimControl$maxit <- 5000
      optRes <- optim(pStartTransform,
                      p_transformed_Laplace3, p_transformed_gr_Laplace3,
                      method = "BFGS", control = optimControl)
      ## Return MLEs on the original scale
      transformedMLEs <- optRes$par
      MLEs <- paramsTransformation$inverseTransform(transformedMLEs)
      optRes$par <- MLEs
      return(optRes)
      returnType(optimResultNimbleList())
   }
  ),
  enableDerivs = list(inverseTransform = list())
)
