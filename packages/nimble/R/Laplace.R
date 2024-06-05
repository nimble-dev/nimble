## NIMBLE Laplace approximation
## AGHQuad/Laplace base class
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
    },
    reset_outer_logLik = function(){},
    save_outer_logLik = function(logLikVal = double()){},
    get_param_value = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
    },
		get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
    },
		get_inner_negHessian = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
    },
		get_inner_negHessian_chol = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
    },
    check_convergence = function(){
      returnType(double())
    },
    updateSettings = function(optimMethod = character(0, default="NULL"),
                               optimStart = character(0, default="NULL"),
                               optimStartValues = double(1, default=Inf),
                               optimWarning = integer(0, default = -1),
                               useInnerCache = integer(0, default=-1),
                               nQuad = integer(0, default=-1),
                               gridType = character(0, default="NULL"),
                               optimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                               replace_optimControl = logical(0, default=FALSE)) {
    },
    ## set_nQuad = function(nQUpdate = integer()){},
    ## set_transformation = function(transformation = character()){},
    ## set_warning = function(warn = logical()){},
    ## set_reInitMethod = function(method = character(), value=double(1)){},
    set_randomeffect_values = function(p = double(1)){}
    ## set_inner_cache = function(cache = logical(0, default = TRUE)){}
  )
)

setup_OneAGHQuad <- function(model, paramNodes, randomEffectsNodes, calcNodes,
                             control) {
  # common setup steps for 1D and >1D cases
  optimControl_ <- extractControlElement(control, 'optimControl', nimOptimDefaultControl())
  optimMethod_ <- extractControlElement(control, 'optimMethod', 'BFGS')
  optimStart_ <- extractControlElement(control, 'optimStart', 'constant')
  optimStartValues_ <- extractControlElement(control, 'optimStartValues', 0)
  nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))

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
  reTrans <- parameterTransform(model, randomEffectsNodes)
  npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))

  if(npar > 1) p_indices <- as.numeric(1:npar)
  else p_indices <- as.numeric(c(1, -1))

  list(optimControl_=optimControl_,
       optimMethod_=optimMethod_,
       optimStart_=optimStart_,
       optimStartValues_=optimStartValues_,
       nre = nre,
       paramDeps = paramDeps,
       innerCalcNodes = innerCalcNodes,
       calcNodes = calcNodes,
       wrtNodes = wrtNodes,
       reTrans = reTrans,
       npar = npar,
       p_indices = p_indices
       )
}

## A single Laplace approximation for only one scalar random effect node
buildOneLaplace1D <- function(model, paramNodes, randomEffectsNodes, calcNodes,
                              control = list()) {
#                              optimControl, optimMethod, optimStart, optimStartValues=0) {
  buildOneAGHQuad1D(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes,
                    control) #optimControl, optimMethod, optimStart, optimStartValues)
}

buildOneAGHQuad1D <- nimbleFunction(
  contains = AGHQuad_BASE,
  setup = function(model, nQuad, paramNodes, randomEffectsNodes, calcNodes,
                   control = list()) {
                   #optimControl, optimMethod, optimStart, optimStartValues=0) {
    ## Check the number of random effects is 1
    ## optimControl_ <- extractControlElement(control, 'optimControl', nimOptimDefaultControl())
    ## optimMethod_ <- extractControlElement(control, 'optimMethod', 'BFGS')
    ## optimStart_ <- extractControlElement(control, 'optimStart', 'constant')
    ## optimStartValues_ <- extractControlElement(control, 'optimStartValues', 0)
    nQuad_ <- nQuad
    S <- setup_OneAGHQuad(model, paramNodes, randomEffectsNodes, calcNodes,
                          control)
    optimControl_ <- S$optimControl_
    optimMethod_ <- S$optimMethod_
    optimStart_ <- S$optimStart_
    optimStartValues_ <- S$optimStartValues_
    nre  <-  S$nre
    paramDeps  <-  S$paramDeps
    innerCalcNodes  <-  S$innerCalcNodes
    calcNodes  <-  S$calcNodes
    wrtNodes  <-  S$wrtNodes
    reTrans  <-  S$reTrans
    npar  <-  S$npar
    p_indices  <-  S$p_indices

    ## nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))
    if(length(nre) != 1) stop("Number of random effects for buildOneAGHQuad1D or buildOneLaplace1D must be 1.")
    ## Check and add necessary upstream deterministic nodes into calcNodes
    ## This ensures that deterministic nodes between paramNodes and calcNodes are used.
    ## paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
    ## if(length(paramDeps) > 0) {
    ##   keep_paramDeps <- logical(length(paramDeps))
    ##   for(i in seq_along(paramDeps)) {
    ##     if(any(paramDeps[i] == calcNodes)) keep_paramDeps[i] <- FALSE
    ##     else {
    ##       nextDeps <- model$getDependencies(paramDeps[i])
    ##       keep_paramDeps[i] <- any(nextDeps %in% calcNodes)
    ##     }
    ##   }
    ##   paramDeps <- paramDeps[keep_paramDeps]
    ## }
    ## innerCalcNodes <- calcNodes
    ## calcNodes <- model$expandNodeNames(c(paramDeps, calcNodes), sort = TRUE)
    ## wrtNodes <- c(paramNodes, randomEffectsNodes)
    ## Indices of randomEffectsNodes and paramNodes inside wrtNodes
    ## npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))
    re_indices <- as.numeric(c(npar+1, -1))
    ## if(npar > 1) p_indices <- as.numeric(1:npar)
    ## else p_indices <- as.numeric(c(1, -1))
    ## ## Indices of randomEffectsNodes inside randomEffectsNodes for use in getting the derivative of
    ## ## the inner log-likelihood (paramNodes fixed) w.r.t. randomEffectsNodes.
    re_indices_inner <- as.numeric(c(1, -1))
    p_and_re_indices <- as.numeric(1:(npar + 1))
    
    ## Set up start values for the inner optimization of Laplace approximation
    if(!is.character(optimStart_) | length(optimStart_) != 1) stop("problem with optimStart ", optimStart_)
    startID <- switch(optimStart_, last=1, last.best=2, constant=3, random=4, model=5)
    if(startID==5) {
      constant_init_par <- c(values(model, randomEffectsNodes), -1)
    } else
      constant_init_par <- c(optimStartValues_, -1)
    
    ## Update and constant nodes for obtaining derivatives using AD
    inner_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes
    
    ## Automated transformation for random effects to ensure range of (-Inf, Inf) 
    ## reTrans <- parameterTransform(model, randomEffectsNodes)
    
    ## The following are used for caching values and gradient in the Laplace3 system
    logLik3_saved_value <- -Inf # numeric(1)
    logLik3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    logLik3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    ## The following are used for caching values for init purposes
    max_inner_logLik_last_argmax <- constant_init_par
    max_inner_logLik_last_value <- -Inf
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    ## Record the maximum Laplace loglikelihood value for obtaining inner optimization start values
    max_logLik <- -Inf
    max_logLik_last_best_argmax <- constant_init_par

    ## Last call cache of neg Hessian.
		saved_inner_negHess <- matrix(0, nrow = 1, ncol = 1)
		## Cache log like saved value to keep track of 3 methods.
    logLik_saved_value <- -Inf
    
    ## Values to save when max inner log lik reached.
    max_outer_logLik <- -Inf
    outer_mode_inner_negHess <- matrix(0, nrow = 1, ncol = 1)
    outer_mode_max_inner_logLik_last_argmax <- as.numeric(c(1, -1))
    outer_param_max <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    
    ## Cached gradients for AGHQ.
    gr_sigmahatwrtre <- numeric(1)
    gr_sigmahatwrtp <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    gr_rehatwrtp <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1)) # double(1)
    gr_QuadSum_value <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    AGHQuad_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(1, -1))
    quadrature_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    
    ## Convergence check for outer function.
    converged <- 0
    
    ## Build AGHQ grid for 1D:
    aghq_grid <- buildAGHQGrid(d = 1, nQuad = nQuad_)
    
    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE    
    
    warn_optim <- extractControlElement(control, 'optimWarning', FALSE) ## Warn about inner optimization issues
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
    one_time_fixes = function() {
      ## Run this once after compiling; remove extraneous -1 if necessary
      if(one_time_fixes_done) return()
      re_indices <<- fix_one_vec(re_indices)
      re_indices_inner <<- fix_one_vec(re_indices_inner)
      max_inner_logLik_last_argmax <<- fix_one_vec(max_inner_logLik_last_argmax)
      outer_mode_max_inner_logLik_last_argmax <<-  fix_one_vec(outer_mode_max_inner_logLik_last_argmax)
      max_logLik_last_best_argmax <<- fix_one_vec(max_logLik_last_best_argmax)
      constant_init_par <<- fix_one_vec(constant_init_par)
      #      if(startID == 3) optStart <<- fix_one_vec(optStart)
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        logLik3_saved_gr <<- fix_one_vec(logLik3_saved_gr)
        logLik3_previous_p <<- fix_one_vec(logLik3_previous_p)
        max_inner_logLik_previous_p <<- fix_one_vec(max_inner_logLik_previous_p)
        outer_param_max <<- fix_one_vec(outer_param_max)
        gr_sigmahatwrtp <<- fix_one_vec(gr_sigmahatwrtp)
        gr_rehatwrtp <<- fix_one_vec(gr_rehatwrtp)
        gr_QuadSum_value <<- fix_one_vec(gr_QuadSum_value)
        AGHQuad_saved_gr <<- fix_one_vec(AGHQuad_saved_gr)
        quadrature_previous_p <<- fix_one_vec(quadrature_previous_p)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit(reInit)
      one_time_fixes_done <<- TRUE
    },
    updateSettings = function(optimMethod = character(0, default="NULL"),
                               optimStart = character(0, default="NULL"),
                               optimStartValues = double(1, default=Inf),
                               optimWarning = integer(0, default = -1),
                               useInnerCache = integer(0, default=-1),
                               nQuad = integer(0, default=-1),
                               gridType = character(0, default="NULL"),
                               optimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                               replace_optimControl = logical(0, default=FALSE)) {
      # Checking should have been done already. Or, if this is being called directly,
      # it will be for development or advanced uses and we can skip checking.
      if(optimMethod != "NULL") optimMethod_ <<- optimMethod
      if(optimStart != "NULL") {
        if(optimStart == "last") startID <<- 1 # last
        else if(optimStart == "last.best") startID <<- 2 # last.best
        else if(optimStart == "constant") startID <<- 3 # use fixed vector optimStart provided at setup time
        else if(optimStart == "random") startID <<- 4
        else if(optimStart == "model") {
          startID <<- 3
          constant_init_par <<- reTrans$transform(values(model, randomEffectsNodes))
        }
      }
      if((length(optimStartValues) != 1) | (optimStartValues[1] != Inf) ) {
        if((length(optimStartValues) == 1) & (optimStartValues[1] == -Inf) ) { # numeric code for "model" setting
          constant_init_par <<- reTrans$transform(values(model, randomEffectsNodes))
        } else {
          if(startID <= 3) {
            constant_init_par <<- optimStartValues
            if(length(constant_init_par) == 1)
              if(nre > 1)
                constant_init_par <<- rep(constant_init_par, nre)
          }
        }
      }
      if(optimWarning != -1) {
        warn_optim <<- optimWarning != 0
      }
      if(useInnerCache != -1) {
        cache_inner_max <<- useInnerCache != 0
      }
      if(nQuad != -1) {
        aghq_grid$setGridSize(nQUpdate = nQuad)
        nQuad_ <<- nQuad
      }
      ## if(gridType != "") {
      ##   transMethod <<- gridType
      ## }
      if(replace_optimControl) {
        optimControl_ <<- optimControl
      }
    },
    set_reInit = function(re = double(1)) {
      reInitTrans <- reTrans$transform(re)
      max_inner_logLik_last_argmax <<- reInitTrans
    },
    get_reInitTrans = function() {
      if(startID == 1) ans <- max_inner_logLik_last_argmax              ## last
      else if(startID == 2) ans <- max_logLik_last_best_argmax    ## last.best
      else if(startID == 3) ans <- constant_init_par                    ## constant
      else if(startID == 4){                                            ## random (prior).
        model$simulate(randomEffectsNodes)
        ans <- reTrans$transform(values(model, randomEffectsNodes))     ## From prior:    
      }
      return(ans)
      returnType(double(1))
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
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik, method = optimMethod_, control = optimControl_)
      if(optRes$convergence != 0 & warn_optim){
        print("Warning: optim did not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      converged <<- optRes$convergence
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Outer check for inner convergence
    check_convergence = function(){
      returnType(double())
      return(converged)
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
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik_internal, method = optimMethod_, control = optimControl_)
      if(optRes$convergence != 0 & warn_optim){
        print("Warning: optim did not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      converged <<- optRes$convergence
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## These two update methods for max_inner_logLik use the same member data caches
    update_max_inner_logLik = function(p = double(1)) {
      optRes <- max_inner_logLik(p)
      max_inner_logLik_last_argmax <<- optRes$par
      max_inner_logLik_last_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_last_argmax)
      returnType(double(1))
    },
    update_max_inner_logLik_internal = function(p = double(1)) {
      optRes <- max_inner_logLik_internal(p)
      max_inner_logLik_last_argmax <<- optRes$par
      max_inner_logLik_last_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_last_argmax)
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
      hess_wrt_p_wrt_re <- joint_logLik_res$hessian[1:npar, npar + 1, 1]      
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
      reTransform <- max_inner_logLik_last_argmax
      maxValue <- max_inner_logLik_last_value
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
      
      if( nQuad_ == 1) {
        ## Laplace Approximation
        logLik_saved_value <<- maxValue - 0.5 * logdetNegHess_value + 0.5 * 1 * log(2*pi)
      }else{
        ## AGHQ Approximation:
        calcLogLik_AGHQ(p)
      }
      logLik3_saved_value <<- logLik_saved_value

      if( nQuad_ == 1 ){
        ## Gradient of Laplace Approx
        AGHQuad_saved_gr <<- gr_logLik_wrt_p - 0.5*(gr_logdetNegHess_wrt_p_v + hess_cross_terms * (gr_logdetNegHess_wrt_re_v / negHessValue))
      }else{
        ## Gradient of AGHQ Approx.
        ## dre_hat/dp = d^2ll/drep / d^2ll/dre^2
        gr_rehatwrtp <<- hess_cross_terms/negHessValue
        ## dsigma_hat/dp (needed at real scale)
        sigma_hat <- 1/sqrt(negHessValue)
        gr_sigmahatwrtp <<- -0.5*gr_logdetNegHess_wrt_p_v*sigma_hat
        gr_sigmahatwrtre <<- -0.5*gr_logdetNegHess_wrt_re_v*sigma_hat
        
        gr_aghq_sum <- gr_AGHQ_nodes(p = p, method = 2) ## Use method 2 for these?
        AGHQuad_saved_gr <<- gr_aghq_sum - 0.5 * (gr_logdetNegHess_wrt_p_v + gr_logdetNegHess_wrt_re_v * gr_rehatwrtp)
      }
      logLik3_saved_gr <<- AGHQuad_saved_gr

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
        max_logLik_last_best_argmax <<- max_inner_logLik_last_argmax
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
    ## Laplace approximation 2: double taping with separate components
    calcLogLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_last_argmax
      maxValue <- max_inner_logLik_last_value
      logdetNegHessian <- logdetNegHess(p, reTransform)
      saved_inner_negHess <<- matrix(exp(logdetNegHessian), nrow = 1, ncol = 1)

      if(nQuad_ == 1){
        ## Laplace approximation.
        logLik_saved_value <<- maxValue - 0.5 * logdetNegHessian + 0.5 * 1 * log(2*pi)
      }else{
        ## Do Quadrature:
        calcLogLik_AGHQ(p)
      }

      if(logLik_saved_value > max_logLik) {
        max_logLik <<- logLik_saved_value
        max_logLik_last_best_argmax <<- max_inner_logLik_last_argmax
      }
      return(logLik_saved_value)
      returnType(double())
    },
    ## Laplace approximation 1: single taping with separate components
    calcLogLik1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_last_argmax
      maxValue <- max_inner_logLik_last_value
      logdetNegHessian <- logdetNegHess(p, reTransform)
      saved_inner_negHess <<- matrix(exp(logdetNegHessian), nrow = 1, ncol = 1)
      
      if(nQuad_ == 1){
        ## Laplace approximation.
        logLik_saved_value <<- maxValue - 0.5 * logdetNegHessian + 0.5 * 1 * log(2*pi)
      }else{
        ## Do Quadrature:
        calcLogLik_AGHQ(p)
      }
      
      if(logLik_saved_value > max_logLik) {
        max_logLik <<- logLik_saved_value
        max_logLik_last_best_argmax <<- max_inner_logLik_last_argmax
      }          
      
      return(logLik_saved_value)
      returnType(double())
    },
    calcLogLik_AGHQ = function(p = double(1)){
      ## AGHQ Approximation:  3 steps. build grid (happens once), transform z to re, save log density.
      aghq_grid$buildGrid()
      nQ <- aghq_grid$getGridSize()
      aghq_grid$transformGrid1D(negHess = saved_inner_negHess, inner_mode = max_inner_logLik_last_argmax)
      modeIndex <- aghq_grid$getModeIndex() ## if even, this is -1
      aghq_grid$saveLogDens( -1, max_inner_logLik_last_value ) ## Cache this value regardless of even or odd.
      for(i in 1:nQ) {
        if(i != modeIndex) aghq_grid$saveLogDens(i, joint_logLik(p = p, reTransform = aghq_grid$getNodesTransformed(i) ) )
      }

      ## Given all the saved values, weights and log density, do quadrature sum.
      logLik_saved_value <<- aghq_grid$quadSum()
      quadrature_previous_p <<- p ## Cache this to make sure you have it for 
    },
    ## Gradient of the Laplace approximation (version 2) w.r.t. parameters
    gr_logLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_last_argmax
      saved_inner_negHess <<- negHess(p, reTransform)
      negHessian <- saved_inner_negHess[1, 1]

      # invNegHessian <- inverse(negHessian)
      grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p(p, reTransform)
      grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re(p, reTransform)[1]
      hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re(p, reTransform)[,1]

      if( nQuad_ == 1 ){
        ## Gradient of Laplace Approx
        p1 <- gr_joint_logLik_wrt_p(p, reTransform)
        AGHQuad_saved_gr <<- p1 - 0.5 * (grlogdetNegHesswrtp + hesslogLikwrtpre * (grlogdetNegHesswrtre / negHessian))
      }else{
        ## Gradient of AGHQ Approx.
        ## dre_hat/dp = d^2ll/drep / d^2ll/dre^2
        gr_rehatwrtp <<- hesslogLikwrtpre/negHessian
        ## dsigma_hat/dp (needed at real scale)
        sigma_hat <- 1/sqrt(negHessian)
        gr_sigmahatwrtp <<- -0.5*grlogdetNegHesswrtp*sigma_hat
        gr_sigmahatwrtre <<- -0.5*grlogdetNegHesswrtre*sigma_hat
        ## Sum gradient of each node.
        gr_aghq_sum <- gr_AGHQ_nodes(p = p, method = 2)
        AGHQuad_saved_gr <<- gr_aghq_sum - 0.5 * (grlogdetNegHesswrtp + grlogdetNegHesswrtre * gr_rehatwrtp)
      }
      return(AGHQuad_saved_gr)
      returnType(double(1))
    },
    ## Gradient of the Laplace approximation (version 1) w.r.t. parameters
    gr_logLik1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_last_argmax
      saved_inner_negHess <<- negHess_internal(p, reTransform)  ## repeated comp. pvdb.
      negHessian <- saved_inner_negHess[1, 1]
      
      ## invNegHessian <- inverse(negHessian)
      grlogdetNegHesswrtp <- gr_logdetNegHess_wrt_p_internal(p, reTransform)
      grlogdetNegHesswrtre <- gr_logdetNegHess_wrt_re_internal(p, reTransform)[1]
      hesslogLikwrtpre <- hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform)[,1]

      if( nQuad_ == 1 ){
        ## Gradient of Laplace Approx
        p1 <- gr_joint_logLik_wrt_p_internal(p, reTransform)
        AGHQuad_saved_gr <<- p1 - 0.5 * (grlogdetNegHesswrtp + hesslogLikwrtpre * (grlogdetNegHesswrtre / negHessian))
      }else{
        ## Gradient of AGHQ Approx.
        ## dre_hat/dp = d^2ll/drep / d^2ll/dre^2
        gr_rehatwrtp <<- hesslogLikwrtpre/negHessian
        ## dsigma_hat/dp (needed at real scale)
        sigma_hat <- 1/sqrt(negHessian)
        gr_sigmahatwrtp <<- -0.5*grlogdetNegHesswrtp*sigma_hat
        gr_sigmahatwrtre <<- -0.5*grlogdetNegHesswrtre*sigma_hat
        ## Sum gradient of each node.
        gr_aghq_sum <- gr_AGHQ_nodes(p = p, method = 1)
        AGHQuad_saved_gr <<- gr_aghq_sum - 0.5 * (grlogdetNegHesswrtp + grlogdetNegHesswrtre * gr_rehatwrtp)
      }
        
      return(AGHQuad_saved_gr)
      returnType(double(1))
    },
    ## Partial gradient of AGHQ nodes w respect to p.
    gr_AGHQ_nodes = function(p = double(1), method = double()){

      ## Need to have quadrature sum for gradient:
      if(any(p != quadrature_previous_p)){
        calcLogLik_AGHQ(p)
      }
 
      ## Method 2 implies double taping.
      modeIndex <- aghq_grid$getModeIndex()
      nQ <- aghq_grid$getGridSize()
      gr_margLogLik_wrt_p <- numeric(value = 0, length = dim(p)[1])
      wgts_lik <- numeric(value = 0, length = nQ)
      for(i in 1:nQ) {
        z_node_i <- aghq_grid$getNodes(i)[1]
        reTrans_i <- aghq_grid$getNodesTransformed(i)
        wgts_lik[i] <- exp(aghq_grid$getLogDensity(i) - max_inner_logLik_last_value)*aghq_grid$getWeights(i)
        
        ## At the mode (z = 0, don't have additional z*sigma_hat gr complication).
	      if( modeIndex == i ){
          if( method == 2 ) gr_jointlogLikwrtp <- gr_joint_logLik_wrt_p(p, reTrans_i)
          else gr_jointlogLikwrtp <- gr_joint_logLik_wrt_p_internal(p, reTrans_i)
          gr_margLogLik_wrt_p <- gr_margLogLik_wrt_p + wgts_lik[i]*gr_jointlogLikwrtp
        }else{
          ## Chain Rule: dll/dre * ( dre_hat/dp + dsigma_hat/dp*z_i )
          ## dll/dp
          if(method == 2){
            gr_logLikwrtrewrtre_i <- gr_joint_logLik_wrt_re(p, reTrans_i)[1]
            gr_logLikewrtp_i <- gr_joint_logLik_wrt_p(p, reTrans_i)
          }else{
            gr_logLikwrtrewrtre_i <- gr_joint_logLik_wrt_re_internal(p, reTrans_i)[1]
            gr_logLikewrtp_i <- gr_joint_logLik_wrt_p_internal(p, reTrans_i)
          }
          gr_logLikwrtrewrtp_i <- gr_logLikwrtrewrtre_i *
                            ( (1 + gr_sigmahatwrtre*z_node_i) * gr_rehatwrtp  +  gr_sigmahatwrtp*z_node_i )
          ## The weighted gradient for the ith sum.
          gr_margLogLik_wrt_p <- gr_margLogLik_wrt_p + wgts_lik[i]*( gr_logLikewrtp_i +  gr_logLikwrtrewrtp_i )
        }
      }
      
      return(gr_margLogLik_wrt_p / sum(wgts_lik[1:nQ]))
      returnType(double(1))
    },
		get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      if(atOuterMode) return(outer_mode_max_inner_logLik_last_argmax)
      return(max_inner_logLik_last_argmax)
    },
		get_inner_negHessian = function(atOuterMode = integer(0, default = 0)){ 
      returnType(double(2))
      if(atOuterMode) return(outer_mode_inner_negHess)
      return(saved_inner_negHess)
    },
		get_inner_negHessian_chol = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
      if(atOuterMode) return(sqrt(outer_mode_inner_negHess))
      return(sqrt(saved_inner_negHess))
    },
    ## Update the maximum mode and neg hess based on the log likelihood passed via optim.
    ##  For efficient saving of values for calculating MLE values of random-effects and INLA simulation of them.
    save_outer_logLik = function(logLikVal = double()){
      if(logLikVal >= max_outer_logLik) {
        max_outer_logLik <<- logLikVal
        outer_mode_inner_negHess <<- saved_inner_negHess
        outer_mode_max_inner_logLik_last_argmax <<- max_inner_logLik_last_argmax
        outer_param_max <<- max_inner_logLik_previous_p
      }
    },
    get_param_value = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      ## Ensures that the inner value will not match and cached values will not be used.
      if(!cache_inner_max) return(numeric(value = Inf, length = npar))
      if(atOuterMode) return(outer_param_max)
      return(max_inner_logLik_previous_p)
    },
    ## Need to reset every time optim is called to recache.
    reset_outer_logLik = function(){
      max_outer_logLik <<- -Inf
    },
    ## Allow the user to explore using different sized quadrature grids.
    ## set_nQuad = function(nQUpdate = integer()){
    ##   aghq_grid$setGridSize(nQUpdate = nQUpdate)
    ##   nQuad <<- nQUpdate
    ## },
    ## set_transformation = function(transformation = character()){}, ## Not applicable to 1 Dimension.
    ## set_warning = function(warn = logical()){
    ##   warn_optim <<- warn
    ## },
    ## Internal option to change initial values.
    ## set_reInitMethod = function(method = character(), values = double(1)) {
    ##   if(method == "last") startID <<- 1 # last
    ##   else if(method == "last.best") startID <<- 2 # last.best
    ##   else if(method == "constant") startID <<- 3 # use fixed vector optimStart provided at setup time
    ##   else if(method == "random") startID <<- 4
    ##   else if(method == "model") {
    ##     startID <<- 3
    ##     constant_init_par <<- reTrans$transform(values(model, randomEffectsNodes))
    ##   } else {
    ##     stop("invalid method for RE initialization")
    ##   }
    ##   if(startID <= 3) {
    ##     constant_init_par <<- values
    ##   }
    ## },
    set_randomeffect_values = function(p = double(1)){
      foundIt <- FALSE
      ## Last value called:
      if(all(p == max_inner_logLik_previous_p)) {
        re <- reTrans$inverseTransform(max_inner_logLik_last_argmax)
        foundIt <- TRUE
      }
      ## Best value called:
      if(all(p == outer_param_max)) {
        re <- reTrans$inverseTransform(outer_mode_max_inner_logLik_last_argmax)
        foundIt <- TRUE
      }
      if(foundIt){
        values(model, paramNodes) <<- p
        model$calculate(paramDeps)
      }else{
        # It would be nice to emit a message here, but different optimizers (e.g. BFGS vs nlminb)
        # behave differently as to whether the previous (last) parameters were always the MLE.
        # print("  [Warning] Have not cached the inner optimization. Running optimization now.")
        update_max_inner_logLik(p)
        re <- reTrans$inverseTransform(max_inner_logLik_last_argmax)
      }
      ## Ensure the model is up to date for all nodes.
      values(model, randomEffectsNodes) <<- re
      model$calculate(innerCalcNodes)
    }
    ## set_inner_cache = function(cache = logical(0, default = TRUE)){
    ##   cache_inner_max <<- cache
    ## }
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
buildOneLaplace <- function(model, paramNodes, randomEffectsNodes, calcNodes,
                            control = list()) {
  #optimControl, optimMethod, optimStart, optimStartValues=0) {
  buildOneAGHQuad(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes,
                  control)
#                  optimControl, optimMethod, optimStart, optimStartValues)
}

buildOneAGHQuad <- nimbleFunction(
  contains = AGHQuad_BASE,
  setup = function(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes,
                   control = list()) {
#                   optimControl, optimMethod, optimStart, optimStartValues=0) {
    ## Check and add necessary (upstream) deterministic nodes into calcNodes
    ## This ensures that deterministic nodes between paramNodes and calcNodes are used.
    ## optimControl_ <- extractControlElement(control, 'optimControl', nimOptimDefaultControl())
    ## optimMethod_ <- extractControlElement(control, 'optimMethod', 'BFGS')
    ## optimStart_ <- extractControlElement(control, 'optimStart', 'constant')
    ## optimStartValues_ <- extractControlElement(control, 'optimStartValues', 0)
    nQuad_ <- nQuad
    S <- setup_OneAGHQuad(model, paramNodes, randomEffectsNodes, calcNodes,
                          control)
    optimControl_ <- S$optimControl_
    optimMethod_ <- S$optimMethod_
    optimStart_ <- S$optimStart_
    optimStartValues_ <- S$optimStartValues_
    nre  <-  S$nre
    paramDeps  <-  S$paramDeps
    innerCalcNodes  <-  S$innerCalcNodes
    calcNodes  <-  S$calcNodes
    wrtNodes  <-  S$wrtNodes
    reTrans  <-  S$reTrans
    npar  <-  S$npar
    p_indices  <-  S$p_indices

    ## paramDeps <- model$getDependencies(paramNodes, determOnly = TRUE, self=FALSE)
    ## if(length(paramDeps) > 0) {
    ##   keep_paramDeps <- logical(length(paramDeps))
    ##   for(i in seq_along(paramDeps)) {
    ##     if(any(paramDeps[i] == calcNodes)) keep_paramDeps[i] <- FALSE
    ##     else {
    ##       nextDeps <- model$getDependencies(paramDeps[i])
    ##       keep_paramDeps[i] <- any(nextDeps %in% calcNodes)
    ##     }
    ##   }
    ##   paramDeps <- paramDeps[keep_paramDeps]
    ## }
    ## innerCalcNodes <- calcNodes
    ## calcNodes <- model$expandNodeNames(c(paramDeps, calcNodes), sort = TRUE)
    ## wrtNodes <- c(paramNodes, randomEffectsNodes)
    ## ## Indices of randomEffectsNodes and paramNodes inside wrtNodes
    ## reTrans <- parameterTransform(model, randomEffectsNodes)
    ## npar <- length(model$expandNodeNames(paramNodes, returnScalarComponents = TRUE))
    ## nre  <- length(model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE))
    nreTrans <- reTrans$getTransformedLength()
    if(nreTrans > 1) reTrans_indices <- as.numeric((npar+1):(npar+nreTrans))
    else reTrans_indices <- as.numeric(c(npar+1, -1)) 
    ## if(npar > 1) p_indices <- as.numeric(1:npar)
    ## else p_indices <- as.numeric(c(1, -1))
    ## ## Indices of randomEffectsNodes inside randomEffectsNodes for use in getting the derivative of
    ## ## the inner log-likelihood (paramNodes fixed) w.r.t. randomEffectsNodes.
    if(nreTrans > 1) reTrans_indices_inner <- as.numeric(1:nreTrans)
    else reTrans_indices_inner <- as.numeric(c(1, -1))
    p_and_reTrans_indices <- as.numeric(1:(npar + nreTrans))
    
    ## Set up start values for the inner optimization of Laplace approximation    
    ## Set up start values for the inner optimization of Laplace approximation
    if(!is.character(optimStart_) | length(optimStart_) != 1) stop("problem with optimStart ", optimStart_)
    startID <- switch(optimStart_, last=1, last.best=2, constant=3, random=4, model=5)
    if(startID==5) {
      constant_init_par <- reTrans$transform(c(values(model, randomEffectsNodes)))
    } else {
      if(length(optimStartValues_) == 1)
        constant_init_par <- rep(optimStartValues_, nreTrans)
      else
        constant_init_par <- optimStartValues_
    }
    if(length(constant_init_par) != nreTrans)
      stop("Wrong length of init values for inner optimization in Laplace or AGHQuad. Have ",
           length(constant_init_par), " but expected ", nreTrans, ".")
    if(length(constant_init_par) == 1) constant_init_par <- c(constant_init_par, -1)

    ## Update and constant nodes info for obtaining derivatives using AD
    inner_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = randomEffectsNodes, calcNodes = innerCalcNodes)
    inner_updateNodes   <- inner_derivsInfo$updateNodes
    inner_constantNodes <- inner_derivsInfo$constantNodes
    joint_derivsInfo    <- makeModelDerivsInfo(model = model, wrtNodes = wrtNodes, calcNodes = calcNodes)
    joint_updateNodes   <- joint_derivsInfo$updateNodes
    joint_constantNodes <- joint_derivsInfo$constantNodes
    
    ## The following are used for caching values and gradient in the Laplace3 system
    logLik3_saved_value <- -Inf #numeric(1)
    logLik3_saved_gr <- if(npar > 1) numeric(npar) else as.numeric(c(0, -1))
    logLik3_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    
    max_inner_logLik_last_argmax <- constant_init_par #if(nreTrans > 1) rep(Inf, nreTrans) else as.numeric(c(Inf, -1))
    max_inner_logLik_last_value <- -Inf #numeric(1)
    max_inner_logLik_previous_p <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))
    cache_inner_max <- TRUE
    
    ## Record the maximum Laplace loglikelihood value for obtaining inner optimization start values
    max_logLik <- -Inf
    max_logLik_last_best_argmax <- constant_init_par #if(nreTrans > 1) rep(Inf, nreTrans) else as.numeric(c(0, -1))
    
    ## The following is used to ensure the one_time_fixes are run when needed.
    one_time_fixes_done <- FALSE
    update_once <- TRUE
    gr_inner_update_once <- TRUE
    gr_inner_logLik_force_update <- TRUE
    gr_inner_logLik_first <- TRUE
    negHess_inner_update_once <- TRUE
    negHess_inner_logLik_force_update <- TRUE
    negHess_inner_logLik_first <- TRUE
    
		## Cache values for access in outer function:
		saved_inner_negHess <- matrix(0, nrow = nre, ncol = nre)
		saved_inner_negHess_chol <- matrix(0, nrow = nre, ncol = nre)
    
    ## Cache log like saved value to keep track of 3 methods.
    logLik_saved_value <- -Inf

    max_outer_logLik <- -Inf
    outer_mode_inner_negHess <- matrix(0, nrow = nre, ncol = nre)
    outer_mode_inner_negHess_chol <- matrix(0, nrow = nre, ncol = nre)
    outer_mode_max_inner_logLik_last_argmax <- if(nreTrans > 1) numeric(nreTrans) else as.numeric(c(0, -1))
    outer_param_max <- if(npar > 1) rep(Inf, npar) else as.numeric(c(Inf, -1))

    ## Build AGHQ grid:
    aghq_grid <- buildAGHQGrid(d = nre, nQuad = nQuad_)
    transMethod <- extractControlElement(control, "gridType", "cholesky")
    
    converged <- 0
    warn_optim <- extractControlElement(control, 'optimWarning', FALSE) ## Warn about inner optimization issues
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
    one_time_fixes = function() {
      if(one_time_fixes_done) return()
      if(nre == 1) {
        reTrans_indices <<- fix_one_vec(reTrans_indices)
        reTrans_indices_inner <<- fix_one_vec(reTrans_indices_inner)
        max_inner_logLik_last_argmax <<- fix_one_vec(max_inner_logLik_last_argmax)
        max_logLik_last_best_argmax <<- fix_one_vec(max_logLik_last_best_argmax)
        constant_init_par <<- fix_one_vec(max_logLik_last_best_argmax)
        outer_mode_max_inner_logLik_last_argmax <<- fix_one_vec(outer_mode_max_inner_logLik_last_argmax)
      }
      if(npar == 1) {
        p_indices <<- fix_one_vec(p_indices)
        logLik3_saved_gr <<- fix_one_vec(logLik3_saved_gr)
        logLik3_previous_p <<- fix_one_vec(logLik3_previous_p)
        max_inner_logLik_previous_p <<- fix_one_vec(max_inner_logLik_previous_p)
        outer_param_max <<- fix_one_vec(outer_param_max)
      }
      reInit <- values(model, randomEffectsNodes)
      set_reInit(reInit)
      one_time_fixes_done <<- TRUE
    },
    updateSettings = function(optimMethod = character(0, default="NULL"),
                               optimStart = character(0, default="NULL"),
                               optimStartValues = double(1, default=Inf),
                               optimWarning = integer(0, default = -1),
                               useInnerCache = integer(0, default=-1),
                               nQuad = integer(0, default=-1),
                               gridType = character(0, default="NULL"),
                               optimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                               replace_optimControl = logical(0, default=FALSE)) {
      # Checking should have been done already. Or, if this is being called directly,
      # it will be for development or advanced uses and we can skip checking.
      if(optimMethod != "NULL") optimMethod_ <<- optimMethod
      if(optimStart != "NULL") {
        if(optimStart == "last") startID <<- 1 # last
        else if(optimStart == "last.best") startID <<- 2 # last.best
        else if(optimStart == "constant") startID <<- 3 # use fixed vector optimStart provided at setup time
        else if(optimStart == "random") startID <<- 4
        else if(optimStart == "model") {
          startID <<- 3
          constant_init_par <<- reTrans$transform(values(model, randomEffectsNodes))
        }
      }
      if((length(optimStartValues) != 1) | (optimStartValues[1] != Inf) ) {
        if((length(optimStartValues) == 1) & (optimStartValues[1] == -Inf) ) { # numeric code for "model" setting
          constant_init_par <<- reTrans$transform(values(model, randomEffectsNodes))
        } else {
          if(startID <= 3) {
            constant_init_par <<- optimStartValues
            if(length(constant_init_par) == 1)
              if(nreTrans > 1)
                constant_init_par <<- rep(constant_init_par, nreTrans)
          }
        }
      }
      if(optimWarning != -1) {
        warn_optim <<- optimWarning != 0
      }
      if(useInnerCache != -1) {
        cache_inner_max <<- useInnerCache != 0
      }
      if(nQuad != -1) {
        aghq_grid$setGridSize(nQUpdate = nQuad)
        nQuad_ <<- nQuad
      }
      if(gridType != "NULL") {
        transMethod <<- gridType
      }
      if(replace_optimControl) {
        optimControl_ <<- optimControl
      }
    },
    set_reInit = function(re = double(1)) {
      reInitTrans <- reTrans$transform(re)
      max_inner_logLik_last_argmax <<- reInitTrans
    },
    get_reInitTrans = function() {
      if(startID == 1) ans <- max_inner_logLik_last_argmax                ## last
      else if(startID == 2) ans <- max_logLik_last_best_argmax            ## last best
      else if(startID == 3) ans <- constant_init_par                  ## constant
      else if(startID == 4){                                            ## random
        model$simulate(randomEffectsNodes)
        ans <- reTrans$transform(values(model, randomEffectsNodes))
      }
      return(ans)
      returnType(double(1))
    },
    ## set_gr_inner_update = function(update = logical(0, default = TRUE)) {
    ##   gr_inner_update_once <<- update
    ## },
    ## set_negHess_inner_update = function(update = logical(0, default = TRUE)) {
    ##   negHess_inner_update_once <<- update
    ## },
    set_params = function(p = double(1)) {
      values(model, paramNodes) <<- p
      model$calculate(paramDeps)
      gr_inner_update_once <<- TRUE
      negHess_inner_update_once <<- TRUE
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
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik, method = optimMethod_, control = optimControl_)
      if(optRes$convergence != 0 & warn_optim){
        print("Warning: optim did not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      converged <<- optRes$convergence
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
      optRes <- optim(reInitTrans, inner_logLik, gr_inner_logLik_internal, method = optimMethod_, control = optimControl_)
      if(optRes$convergence != 0 & warn_optim){
        print("Warning: optim did not converge for the inner optimization of AGHQuad or Laplace approximation")
      }
      converged <<- optRes$convergence
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## Outer check on innner convergence.
    check_convergence = function(){
      returnType(double())
      return(converged)
    },    
    ## These two update methods for max_inner_logLik use the same member data caches
    update_max_inner_logLik = function(p = double(1)) {
      optRes <- max_inner_logLik(p)
      max_inner_logLik_last_argmax <<- optRes$par
      max_inner_logLik_last_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_last_argmax)
      returnType(double(1))
    },
    update_max_inner_logLik_internal = function(p = double(1)) {
      optRes <- max_inner_logLik_internal(p)
      max_inner_logLik_last_argmax <<- optRes$par
      max_inner_logLik_last_value <<- optRes$value
      max_inner_logLik_previous_p <<- p
      return(max_inner_logLik_last_argmax)
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
    cholNegHessian = function(p = double(1), reTransform = double(1)) {
      negHessian <- negHess(p, reTransform)
      ans <- chol(negHessian)
      return(ans)
      returnType(double(2))
    },
    ## Logdet negative Hessian
    logdetNegHess = function(p = double(1), reTransform = double(1)) {
      ans <- 2 * sum(log(diag(cholNegHessian(p, reTransform))))
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
      # for(i in 1:nreTrans) negHessUpper[i,i:nreTrans] <- -joint_logLik_res$hessian[npar + i, npar + i:nreTrans, 1]
      cholNegHess <- chol(negHessUpper)
      logdetNegHessAns <- 2 * sum(log(diag(cholNegHess)))
      hess_wrt_p_wrt_re <- matrix(init = FALSE, nrow = npar, ncol = nre)
      for(i in 1:npar){
        for(j in 1:nreTrans){
          hess_wrt_p_wrt_re[i, j] <- joint_logLik_res$hessian[i, npar + j, 1]
        }
      }
      # hess_wrt_p_wrt_re <- joint_logLik_res$hessian[1:npar, npar + (1:nreTrans), 1] # Wasn't working.
      
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
      reTransform <- max_inner_logLik_last_argmax
      maxValue <- max_inner_logLik_last_value
      ans <- derivs(joint_logLik_with_higher_derivs(p, reTransform), wrt = p_and_reTrans_indices, order = 0, model = model,
                    updateNodes = joint_updateNodes, constantNodes = joint_constantNodes)
      ind <- 1
      # all "logLik" here is joint log likelihood (i.e. for p and re)
      gr_logLik_wrt_p <- ans$value[(ind):(ind + npar - 1)]
      ind <- ind + npar
      logdetNegHess_value <- ans$value[ind]
      ind <- ind + 1
      chol_negHess <- matrix(ans$value[(ind):(ind + nreTrans*nreTrans - 1)], nrow = nreTrans, ncol = nreTrans)
      saved_inner_negHess_chol <<- chol_negHess ## Method 3 doesn't cache neg Hessian.*** Should we calc here?
      ind <- ind + nreTrans*nreTrans
      hess_cross_terms <- matrix(ans$value[(ind):(ind + npar*nreTrans - 1)], nrow = npar, ncol = nreTrans)
      ind <- ind + npar*nreTrans
      gr_logdetNegHess_wrt_p_v <- ans$value[(ind):(ind + npar - 1)]
      ind <- ind + npar
      gr_logdetNegHess_wrt_re_v <- ans$value[(ind):(ind + nreTrans - 1)]
      
      if( nQuad_ == 1) {
        ## Laplace Approximation
        logLik_saved_value <<- maxValue - 0.5 * logdetNegHess_value + 0.5 * nreTrans * log(2*pi)
      }else{
        ## AGHQ Approximation:
        calcLogLik_AGHQ(p)
      }
      logLik3_saved_value <<- logLik_saved_value
            
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
        max_logLik_last_best_argmax <<- max_inner_logLik_last_argmax
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
    ## Laplace approximation 2: double taping with separate components
    calcLogLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_last_argmax
      maxValue <- max_inner_logLik_last_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed
      saved_inner_negHess_chol <<- cholNegHessian(p, reTransform)
      logdetNegHessian <- 2 * sum(log(diag(saved_inner_negHess_chol)))

      if(nQuad_ == 1){
        logLik_saved_value <<- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      }else{
        calcLogLik_AGHQ(p)
      }
      if(logLik_saved_value > max_logLik) {
        max_logLik <<- logLik_saved_value
        max_logLik_last_best_argmax <<- max_inner_logLik_last_argmax
      }

      return(logLik_saved_value)
      returnType(double())
    },
    ## Laplace approximation 1: single taping with separate components
    calcLogLik1 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik_internal(p)
      }
      reTransform <- max_inner_logLik_last_argmax
      maxValue <- max_inner_logLik_last_value
      if(maxValue == -Inf) return(-Inf) # This would mean inner optimization failed
      saved_inner_negHess_chol <<- cholNegHessian(p, reTransform)
      logdetNegHessian <- 2 * sum(log(diag(saved_inner_negHess_chol)))
      
      if(nQuad_ == 1){
        ## Laplace Approx
        logLik_saved_value <<- maxValue - 0.5 * logdetNegHessian + 0.5 * nreTrans * log(2*pi)
      }else{
        ## AGHQ Approx
        calcLogLik_AGHQ(p)
      }
      if(logLik_saved_value > max_logLik) {
        max_logLik <<- logLik_saved_value
        max_logLik_last_best_argmax <<- max_inner_logLik_last_argmax
      }
      return(logLik_saved_value)
      returnType(double())
    },
    calcLogLik_AGHQ = function(p = double(1)){
      ## AGHQ Approximation:  3 steps. build grid (happens once), transform z to re, save log density.
      aghq_grid$buildGrid()
      aghq_grid$transformGrid(cholNegHess = saved_inner_negHess_chol, 
                                    inner_mode = max_inner_logLik_last_argmax, method = transMethod)
      modeIndex <- aghq_grid$getModeIndex()
      nQ <- aghq_grid$getGridSize()
      aghq_grid$saveLogDens(-1, max_inner_logLik_last_value )
      for(i in 1:nQ) {
        if(i != modeIndex) aghq_grid$saveLogDens(i, joint_logLik(p = p, reTransform = aghq_grid$getNodesTransformed(i) ) )
      }
      ## Given all the saved values, weights and log density, do quadrature sum.
      logLik_saved_value <<- aghq_grid$quadSum()
    },    
    ## Gradient of the Laplace approximation 2 w.r.t. parameters
    gr_logLik2 = function(p = double(1)){
      if(!one_time_fixes_done) one_time_fixes()
      if(any(p != max_inner_logLik_previous_p) | !cache_inner_max) {
        update_max_inner_logLik(p)
      }
      reTransform <- max_inner_logLik_last_argmax
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
      reTransform <- max_inner_logLik_last_argmax
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
    get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      if(atOuterMode) return(outer_mode_max_inner_logLik_last_argmax)
      return(max_inner_logLik_last_argmax)
    },
		get_inner_negHessian = function(atOuterMode = integer(0, default = 0)){ 
      returnType(double(2))
      if(atOuterMode) return(outer_mode_inner_negHess)
      return(saved_inner_negHess)
    },
		get_inner_negHessian_chol = function(atOuterMode = integer(0, default = 0)){
      returnType(double(2))
      if(atOuterMode) return(outer_mode_inner_negHess_chol)
      return(saved_inner_negHess_chol)
    },
    ## Update the maximum mode and neg hess based on the log likelihood passed via optim.
    ## For efficient saving of values for calculating MLE values of random-effects.
    save_outer_logLik = function(logLikVal = double()){
      if(logLikVal > max_outer_logLik) {
        max_outer_logLik <<- logLikVal
        outer_mode_inner_negHess <<- saved_inner_negHess
        outer_mode_max_inner_logLik_last_argmax <<- max_inner_logLik_last_argmax
        outer_mode_inner_negHess_chol <<- saved_inner_negHess_chol
        outer_param_max <<- max_inner_logLik_previous_p
      }
    },
    get_param_value = function(atOuterMode = integer(0, default = 0)){
      returnType(double(1))
      ## Ensures that the inner value will not match and cached values will not be used.
      if(!cache_inner_max) return(numeric(value = Inf, length = npar))
      if(atOuterMode) return(outer_param_max)
      return(max_inner_logLik_previous_p)
    },    
    ## Need to reset every call optim to recache.
    reset_outer_logLik = function(){
      max_outer_logLik <<- -Inf
    },
    ## set_nQuad = function(nQUpdate = integer()){
    ##   aghq_grid$setGridSize(nQUpdate = nQUpdate)
    ##   nQuad <<- nQUpdate
    ## },
    ## Choose spectral vs cholesky.
    ## set_transformation = function(transformation = character()){
    ##   transMethod <<- transformation
    ## },
    ## set_warning = function(warn = logical()){
    ##   warn_optim <<- warn
    ## },
    ## set_reInitMethod = function(method = character(), values = double(1)) {
    ##   if(method == "last") startID <<- 1 # last
    ##   else if(method == "last.best") startID <<- 2 # last.best
    ##   else if(method == "constant") startID <<- 3 # use fixed vector optimStart provided at setup time
    ##   else if(method == "random") startID <<- 4
    ##   else if(method == "model") {
    ##     startID <<- 3
    ##     constant_init_par <<- reTrans$transform(values(model, randomEffectsNodes))
    ##   } else {
    ##     stop("invalid method for RE initialization")
    ##   }
    ##   if(startID <= 3) {
    ##     constant_init_par <<- values
    ##     if(length(values) == 1)
    ##       if(nreTrans > 1)
    ##         constant_init_par <<- rep(values, nreTrans)
    ##   }
    ## },
    set_randomeffect_values = function(p = double(1)){
      foundIt <- FALSE
      ## Last value called:
      if(all(p == max_inner_logLik_previous_p)) {
        re <- reTrans$inverseTransform(max_inner_logLik_last_argmax)
        foundIt <- TRUE
      }
      ## Best value called:
      if(all(p == outer_param_max)) {
        re <- reTrans$inverseTransform(outer_mode_max_inner_logLik_last_argmax)
        foundIt <- TRUE
      }
      if(foundIt){
        values(model, paramNodes) <<- p
        ans <- model$calculate(paramDeps)
      }else{
        # It would be nice to emit a message here, but different optimizers (e.g. BFGS vs nlminb)
        # behave differently as to whether the previous (last) parameters were always the MLE.
        # print("  [Warning] Have not cached the inner optimization. Running optimization now.")
        update_max_inner_logLik(p)
        re <- reTrans$inverseTransform(max_inner_logLik_last_argmax)
      }
      ## Ensure the model is up to date for all nodes.
      values(model, randomEffectsNodes) <<- re
      model$calculate(innerCalcNodes)
    }
    ## set_inner_cache = function(cache = logical(0, default = TRUE)){
    ##   cache_inner_max <<- cache
    ## }
  ),
  buildDerivs = list(inner_logLik                            = list(),
                     joint_logLik                            = list(),
                     gr_joint_logLik_wrt_re                  = list(),
                     negHess                                 = list(),
                     cholNegHessian                          = list(),
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

#' Organize model nodes for marginalization
#'
#' Process model to organize nodes for marginalization (integration over latent 
#' nodes or random effects) as by Laplace approximation.
#'
#' @param model A nimble model such as returned by \code{nimbleModel}.
#'
#' @param paramNodes A character vector of names of stochastic nodes that are
#'   parameters of nodes to be marginalized over (\code{randomEffectsNodes}).
#'   See details for default.
#'
#' @param randomEffectsNodes A character vector of nodes to be marginalized over
#'   (or "integrated out"). In the case of calculating the likelihood of a model
#'   with continuous random effects, the nodes to be marginalized over are the
#'   random effects, hence the name of this argument. However, one can
#'   marginalize over any nodes desired as long as they are continuous. 
#'   See details for default.
#'
#' @param calcNodes A character vector of nodes to be calculated as the
#'   integrand for marginalization. Typically this will include
#'   \code{randomEffectsNodes} and some data nodes. Se details for default.
#'
#' @param calcNodesOther A character vector of nodes to be calculated as part of
#'   the log likelihood that are not connected to the \code{randomEffectNodes}
#'   and so are not actually part of the marginalization. These are somewhat
#'   extraneous to the purpose of this function, but it is convenient to handle
#'   them here because often the purpose of marginalization is to calculate log
#'   likelihoods, including from "other" parts of the model.
#'
#' @param split A logical indicating whether to split \code{randomEffectsNodes}
#'   into conditionally independent sets that can be marginalized separately
#'   (\code{TRUE}) or to keep them all in one set for a single marginalization
#'   calculation.
#'
#' @param check A logical indicating whether to try to give reasonable warnings
#'   of badly formed inputs that might be missing important nodes or include
#'   unnecessary nodes.
#'
#' @param allowDiscreteLatent A logical indicating whether to
#'   allow discrete latent states. (default = \code{FALSE})
#'
#' @details
#'
#' This function is used by \code{buildLaplace} to organize model nodes into
#' roles needed for setting up the (approximate) marginalization done by Laplace
#' approximation. It is also possible to call this function directly and pass
#' the resulting list (possibly modified for your needs) to \code{buildLaplace}.
#'
#' Any of the input node vectors, when provided, will be processed using
#'   \code{nodes <- model$expandNodeNames(nodes)}, where \code{nodes} may be
#'   \code{paramNodes}, \code{randomEffectsNodes}, and so on. This step allows
#'   any of the inputs to include node-name-like syntax that might contain
#'   multiple nodes. For example, \code{paramNodes = 'beta[1:10]'} can be
#'   provided if there are actually 10 scalar parameters, 'beta[1]' through
#'   'beta[10]'. The actual node names in the model will be determined by the
#'   \code{exapndNodeNames} step.
#'
#' This function does not do any of the marginalization calculations. It only
#' organizes nodes into roles of parameters, random effects, integrand
#' calculations, and other log likelihood calculations.
#'
#' The checking done if `check=TRUE` tries to be reasonable, but it can't cover
#' all cases perfectly. If it gives an unnecessary warning, simply set `check=FALSE`.
#'
#' If \code{paramNodes} is not provided, its default depends on what other
#'   arguments were provided. If neither \code{randomEffectsNodes} nor
#'   \code{calcNodes} were provided, \code{paramNodes} defaults to all
#'   top-level, stochastic nodes, excluding any posterior predictive nodes
#'   (those with no data anywhere downstream). These are determined by
#'   \code{model$getNodeNames(topOnly = TRUE, stochOnly = TRUE,
#'   includePredictive = FALSE)}. If \code{randomEffectsNodes} was provided,
#'   \code{paramNodes} defaults to stochastic parents of
#'   \code{randomEffectsNodes}. In these cases, any provided \code{calcNodes} or
#'   \code{calcNodesOther} are excluded from default \code{paramNodes}. If
#'   \code{calcNodes} but not \code{randomEffectsNodes} was provided, then the
#'   default for \code{randomEffectsNodes} is determined first, and then
#'   \code{paramNodes} defaults to stochastic parents of
#'   \code{randomEffectsNodes}. Finally, any stochastic parents of
#'   \code{calcNodes} (whether provided or default) that are not in
#'   \code{calcNodes} are added to the default for \code{paramNodes}, but only
#'   after \code{paramNodes} has been used to determine the defaults for
#'   \code{randomEffectsNodes}, if necessary.
#'
#' Note that to obtain sensible defaults, some nodes must have been marked as
#'   data, either by the \code{data} argument in \code{nimbleModel} or by
#'   \code{model$setData}. Otherwise, all nodes will appear to be posterior
#'   predictive nodes, and the default \code{paramNodes} may be empty.
#'
#' For purposes of \code{buildLaplace}, \code{paramNodes} does not need to (but
#'   may) include deterministic nodes between the parameters and any
#'   \code{calcNodes}. Such deterministic nodes will be included in
#'   calculations automatically when needed.
#'
#' If \code{randomEffectsNodes} is missing, the default is a bit complicated: it
#'   includes all latent nodes that are descendants (or "downstream") of
#'   \code{paramNodes} (if provided) and are either (i) ancestors (or
#'   "upstream") of data nodes (if \code{calcNodes} is missing), or (ii)
#'   ancestors or elements of \code{calcNodes} (if \code{calcNodes} and
#'   \code{paramNodes} are provided), or (iii) elements of \code{calcNodes} (if
#'   \code{calcNodes} is provided but \code{paramNodes} is missing). In all
#'   cases, discrete nodes (with warning if \code{check=TRUE}), posterior
#'   predictive nodes and \code{paramNodes} are excluded.
#'
#' \code{randomEffectsNodes} should only include stochastic nodes.
#'
#' If \code{calcNodes} is missing, the default is \code{randomEffectsNodes} and
#'   their descendants to the next stochastic nodes, excluding posterior
#'   predictive nodes. These are determined by
#'   \code{model$getDependencies(randomEffectsNodes, includePredictive=FALSE)}.
#'
#' If \code{calcNodesOther} is missing, the default is all stochastic
#'   descendants of \code{paramNodes}, excluding posterior predictive nodes
#'   (from \code{model$getDependencies(paramNodes, stochOnly=TRUE, self=FALSE,
#'   includePosterior=FALSE)}) that are not part of \code{calcNodes}.
#'
#' For purposes of \code{buildLaplace}, neither \code{calcNodes} nor
#'   \code{calcNodesOther} needs to (but may) contain deterministic nodes
#'   between \code{paramNodes} and \code{calcNodes} or \code{calcNodesOther},
#'   respectively. These will be included in calculations automatically when
#'   needed.
#'
#' If \code{split} is \code{TRUE}, \code{model$getConditionallyIndependentSets}
#'   is used to determine sets of the \code{randomEffectsNodes} that can be
#'   independently marginalized. The \code{givenNodes} are the
#'   \code{paramNodes} and \code{calcNodes} excluding any
#'   \code{randomEffectsNodes} and their deterministic descendants. The
#'   \code{nodes} (to be split into sets) are the \code{randomEffectsNodes}.
#'
#' If \code{split} is a numeric vector, \code{randomEffectsNodes} will be split
#'   by \code{split}(\code{randomEffectsNodes}, \code{control$split}). The last
#'   option allows arbitrary control over how \code{randomEffectsNodes} are
#'   blocked.
#'
#' If \code{check=TRUE}, then defaults for each of the four categories of nodes
#'   are created even if the corresponding argument was provided. Then warnings
#'   are emitted if there are any extra (potentially unnecessary) nodes provided
#'   compared to the default or if there are any nodes in the default that were
#'   not provided (potentially necessary). These checks are not perfect and may
#'   be simply turned off if you are confident in your inputs.
#'
#' (If \code{randomEffectsNodes} was provided but \code{calcNodes} was not
#'   provided, the default (for purposes of \code{check=TRUE} only) for
#'   \code{randomEffectsNodes} differs from the above description. It uses
#'   stochastic descendants of \code{randomEffectsNodes} in place of the
#'   "data nodes" when determining ancestors of data nodes. And it uses item
#'   (ii) instead of (iii) in the list above.)
#'
#' @author Wei Zhang, Perry de Valpine, Paul van Dam-Bates
#' @return
#'
#' A list is returned with elements:
#'
#' \itemize{
#'
#' \item \code{paramNodes}: final processed version of \code{paramNodes}
#'
#' \item \code{randomEffectsNodes}: final processed version of \code{randomEffectsNodes}
#'
#' \item \code{calcNodes}: final processed version of \code{calcNodes}
#'
#' \item \code{calcNodesOther}: final processed version of \code{calcNodesOther}
#'
#' \item \code{givenNodes}: Input to \code{model$getConditionallyIndependentSets}, if \code{split=TRUE}.
#'
#' \item \code{randomEffectsSets}: Output from
#'   \code{model$getConditionallyIndependentSets}, if \code{split=TRUE}. This
#'   will be a list of vectors of node names. The node names in one list element
#'   can be marginalized independently from those in other list elements. The
#'   union of the list elements should be all of \code{randomEffectsNodes}. If
#'   \code{split=FALSE}, \code{randomEffectsSets} will be a list with one
#'   element, simply containing \code{randomEffectsNodes}. If \code{split} is a
#'   numeric vector,  \code{randomEffectsSets} will be the result of
#'   \code{split}(\code{randomEffectsNodes}, \code{control$split}).
#'
#' }
#'
#' @export
setupMargNodes <- function(model, paramNodes, randomEffectsNodes, calcNodes,
                           calcNodesOther,
                           split = TRUE,
                           check = TRUE,
                           allowDiscreteLatent = FALSE) {
  paramProvided     <- !missing(paramNodes)
  reProvided        <- !missing(randomEffectsNodes)
  calcProvided      <- !missing(calcNodes)
  calcOtherProvided <- !missing(calcNodesOther)

  normalizeNodes <- function(nodes, sort = FALSE) {
    if(is.null(nodes) || isFALSE(nodes)) character(0)
    else model$expandNodeNames(nodes, sort = sort)
  }
  if(paramProvided) paramNodes         <- normalizeNodes(paramNodes)
  if(reProvided)    randomEffectsNodes <- normalizeNodes(randomEffectsNodes)
  if(calcProvided)  calcNodes          <- normalizeNodes(calcNodes, sort = TRUE)
  if(calcOtherProvided) calcNodesOther <- normalizeNodes(calcNodesOther, sort = TRUE)

  if(reProvided) {
    if(check && !allowDiscreteLatent)
      if(any(model$isDiscrete(randomEffectsNodes)))
        warning("Some randomEffectsNodes follow discrete distributions. That is likely to cause problems.")
  }

  # We considered a feature to allow params to be nodes without priors. This is a placeholder in case
  # we ever pursue that again.
  # allowNonPriors <- FALSE
  # We may need to use determ and stochastic dependencies of parameters multiple times below
  # Define these to avoid repeated computation
  # A note for future: determ nodes between parameters and calcNodes are needed inside buildOneAGHQuad
  # and buildOneAGHQuad1D. In the future, these could be all done here to be more efficient
  paramDetermDeps <- character(0)
  paramStochDeps  <- character(0)
  paramDetermDepsCalculated <- FALSE
  paramStochDepsCalculated  <- FALSE

  # 1. Default parameters are stochastic top-level nodes. (We previously
  #    considered an argument allowNonPriors, defaulting to FALSE. If TRUE, the
  #    default params would be all top-level stochastic nodes with no RHSonly
  #    nodes as parents and RHSonly nodes (handling of constants TBD, since
  #    non-scalars would be converted to data) that have stochastic dependencies
  #    (And then top-level stochastic nodes with RHSonly nodes as parents are
  #    essentially latent/data nodes, some of which would need to be added to
  #    randomEffectsNodes below.) However this got too complicated. It is
  #    simpler and clearer to require "priors" for parameters, even though prior
  #    probs may not be used.
  paramsHandled <- TRUE
  if(!paramProvided) {
    if(!reProvided) {
      if(!calcProvided) {
        paramNodes <- model$getNodeNames(topOnly = TRUE, stochOnly = TRUE, includePredictive = FALSE)
      } else {
        # calcNodes were provided, but RE nodes were not, so delay creating default params
        paramsHandled <- FALSE
      }
    } else {
      nodesToFindParentsFrom <- randomEffectsNodes
      paramNodes <- model$getParents(nodesToFindParentsFrom, self=FALSE, stochOnly=TRUE)
      # self=FALSE doesn't omit if one RE node is a parent of another, so we have to do the next step
      paramNodes <- setdiff(paramNodes, nodesToFindParentsFrom)
    }
    if(paramsHandled) {
      if(calcProvided) paramNodes <- setdiff(paramNodes, calcNodes)
      if(calcOtherProvided) paramNodes <- setdiff(paramNodes, calcNodesOther)
    }
  }

  # 2. Default random effects are latent nodes that are downstream stochastic dependencies of params.
  #    In step 3, default random effects are also limited to those that are upstream parents of calcNodes
  if((!reProvided) || check) {
    latentNodes <- model$getNodeNames(latentOnly = TRUE, stochOnly = TRUE,
                                      includeData = FALSE, includePredictive = FALSE)
    if(!allowDiscreteLatent) {
      latentDiscrete <- model$isDiscrete(latentNodes)
      if(any(latentDiscrete)) {
        if((!reProvided) && check) {
          warning("In trying to determine default randomEffectsNodes, there are some nodes\n",
                  "that follow discrete distributions. These will be omitted.")
        }
        latentNodes <- latentNodes[!latentDiscrete]
      }
    }
    if(paramsHandled) {
      paramDownstream <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE,
                                               downstream = TRUE, includePredictive = FALSE)
      #    paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, self = FALSE)
      #    paramStochDepsCalculated <- TRUE
      reNodesDefault <- intersect(latentNodes, paramDownstream)
    } else {
      reNodesDefault <- latentNodes
    }
    # Next, if calcNodes were not provided, we create a temporary
    # dataNodesDefault for purposes of updating reNodesDefault if needed. The
    # idea is that reNodesDefault should be trimmed to include only nodes
    # upstream of "data" nodes, where "data" means nodes in the role of data for
    # purposes of marginalization.
    # The tempDataNodesDefault is either dependencies of RE nodes if provided, or
    # actual data nodes in the model if RE nodes not provided.
    # If calcNodes were provided, then they are used directly to trim reNodesDefault.
    if(!calcProvided) {
      if(reProvided)
        tempDataNodesDefault <- model$getDependencies(randomEffectsNodes, stochOnly = TRUE,
                                                      self = FALSE, includePredictive = FALSE)
      else
        tempDataNodesDefault <- model$getNodeNames(dataOnly = TRUE)
      if(paramsHandled)
        tempDataNodesDefault <- setdiff(tempDataNodesDefault, paramNodes)
      tempDataNodesDefaultParents <- model$getParents(tempDataNodesDefault, upstream = TRUE, stochOnly = TRUE)
      # See comment above about why this is necessary:
      tempDataNodesDefaultParents <- setdiff(tempDataNodesDefaultParents, tempDataNodesDefault)
      reNodesDefault <- intersect(reNodesDefault, tempDataNodesDefaultParents)
    } else {
      # Update reNodesDefault to exclude nodes that lack downstream connection to a calcNode
      if(paramsHandled) { # This means reProvided OR paramsProvided. Including parents allows checking
        # of potentially missing REs.
        reNodesDefault <- intersect(reNodesDefault,
                                    model$getParents(calcNodes, upstream=TRUE, stochOnly = TRUE))
      } else { # This means !paramsHandled and hence !reProvided AND !paramsProvided
        reNodesDefault <- intersect(reNodesDefault,
                                    calcNodes)
        reNodesDefault <- intersect(reNodesDefault,
                                    model$getParents(calcNodes, upstream=TRUE, stochOnly = TRUE))
      }
    }
  }

  # If only calcNodes were provided, we have now created reNodesDefault from calcNodes,
  # and are now ready to create default paramNodes
  if(!paramsHandled) {
    paramNodes <- model$getParents(reNodesDefault, self=FALSE, stochOnly=TRUE)
    # See comment above about why this is necessary:
    paramNodes <- setdiff(paramNodes, reNodesDefault)
    if(calcOtherProvided) paramNodes <- setdiff(paramNodes, calcNodesOther)
  }

  # 3. Optionally check random effects if they were provided (not default)
  if(reProvided && check) {
    # First check is for random effects that should have been included but weren't
    reCheck <- setdiff(reNodesDefault, randomEffectsNodes)
    if(length(reCheck)) {
      errorNodes <- paste0(head(reCheck, n = 4), sep = "", collapse = ", ")
      if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
      warning(paste0("There are some random effects (latent states) in the model that look\n",
                     "like they should be included in randomEffectsNodes for Laplace or AGHQuad approximation\n",
                     "for the provided (or default) paramNodes:\n",
                     errorNodes, "\n",
                     "To silence this warning, include \'check = FALSE\' in the control list\n",
                     "to buildLaplace or as an argument to setupMargNodes."))
    }
    # Second check is for random effects that were included but look unnecessary
    reCheck <- setdiff(randomEffectsNodes, reNodesDefault)
    if(length(reCheck)) {
      # Top nodes should never trigger warning.
      # Descendants of top nodes that are in randomEffectsNodes should not trigger warning
      topNodes <- model$getNodeNames(topOnly=TRUE)
      reCheckTopNodes <- intersect(reCheck, topNodes)
      if(length(reCheckTopNodes)) {
        # Simple downstream=TRUE here is not a perfect check of connection among all nodes
        # but it will avoid false alarms
        reCheck <- setdiff(reCheck, model$getDependencies(reCheckTopNodes, downstream=TRUE, stochOnly=TRUE))
      }
      if(length(reCheck)) {
        errorNodes <- paste0(head(reCheck, n = 4), sep = "", collapse = ", ")
        if(length(reCheck) > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some randomEffectsNodes provided that look like\n",
                       "they are not needed for Laplace or AGHQuad approximation for the\n",
                       "provided (or default) paramNodes:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'check = FALSE\' in the control list\n",
                       "to buildLaplace or as an argument to setupMargNodes."))
      }
    }
  }
  # Set final choice of randomEffectsNodes
  if(!reProvided) {
    randomEffectsNodes <- reNodesDefault
  }

  # Set actual default calcNodes. This time it has self=TRUE (default)
  if((!calcProvided) || check) {
    calcNodesDefault <- model$getDependencies(randomEffectsNodes, includePredictive = FALSE)
  }
  # 5. Optionally check calcNodes if they were provided (not default)
  if(calcProvided && check) {
    # First check is for calcNodes that look necessary but were omitted
    calcCheck <- setdiff(calcNodesDefault, calcNodes)
    if(length(calcCheck)) {
      errorNodes <- paste0(head(calcCheck, n = 4), sep = "", collapse = ", ")
      if(length(calcCheck) > 4) errorNodes <- paste(errorNodes, "...")
      warning(paste0("There are some model nodes that look like they should be\n",
                     "included in the calcNodes for Laplace or AGHQuad approximation because\n",
                     "they are dependencies of some randomEffectsNodes:\n",
                     errorNodes, "\n",
                     "To silence this warning, include \'check = FALSE\' in the control list\n",
                     "to buildLaplace or as an argument to setupMargNodes."))
    }
    # Second check is for calcNodes that look unnecessary
    # If some determ nodes between paramNodes and randomEffectsNodes are provided in calcNodes
    # then that's ok and we should not throw a warning message.
    calcCheck <- setdiff(calcNodes, calcNodesDefault)
    errorNodes <- calcCheck[model$getNodeType(calcCheck)=="stoch"]
    # N.B. I commented out this checking of deterministic nodes for now.
    #      Iterating through individual nodes for getDependencies can be slow
    #      and I'd like to think more about how to do this. -Perry
    ## determCalcCheck <- setdiff(calcCheck, errorNodes)
    ## lengthDetermCalcCheck <- length(determCalcCheck)
    ## # Check other determ nodes
    ## if(lengthDetermCalcCheck){
    ##   paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE, includePredictive = FALSE)
    ##   paramDetermDepsCalculated <- TRUE
    ##   for(i in 1:lengthDetermCalcCheck){
    ##     if(!(determCalcCheck[i] %in% paramDetermDeps) ||
    ##        !(any(model$getDependencies(determCalcCheck[i], self = FALSE) %in% calcNodesDefault))){
    ##       errorNodes <- c(errorNodes, determCalcCheck[i])
    ##     }
    ##   }
    ## }
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      warning(paste0("There are some calcNodes provided that look like\n",
                     "they are not needed for Laplace or AGHQuad approximation over\n",
                     "the provided (or default) randomEffectsNodes:\n",
                     outErrorNodes, "\n",
                     "To silence this warning, include \'check = FALSE\' in the control list\n",
                     "to buildLaplace or as an argument to setupMargNodes."))
    }
  }
  # Finish step 4
  if(!calcProvided){
    calcNodes <- calcNodesDefault
  }
  if(!paramProvided) {
    possibleNewParamNodes <- model$getParents(calcNodes, self=FALSE, stochOnly=TRUE)
    # self=FALSE doesn't omit if one node is a parent of another, so we have to do the next step
    possibleNewParamNodes <- setdiff(possibleNewParamNodes, calcNodesDefault)
    paramNodes <- unique(c(paramNodes, possibleNewParamNodes))
  }

  # 6. Default calcNodesOther: nodes needed for full model likelihood but
  #    that are not involved in the marginalization done by Laplace.
  #    Default is a bit complicated: All dependencies from paramNodes to
  #    stochastic nodes that are not part of calcNodes. Note that calcNodes
  #    does not necessarily contain deterministic nodes between paramNodes and
  #    randomEffectsNodes. We don't want to include those in calcNodesOther.
  #    (A deterministic that is needed for both calcNodes and calcNodesOther should be included.)
  #    So we have to first do a setdiff on stochastic nodes and then fill in the
  #    deterministics that are needed.
  if(!calcOtherProvided || check) {
    paramStochDeps <- model$getDependencies(paramNodes, stochOnly = TRUE, # Should this be dataOnly=TRUE?
                                            self = FALSE, includePredictive = FALSE)
    calcNodesOtherDefault <- setdiff(paramStochDeps, calcNodes)
  }
  if(calcOtherProvided) {
    if((length(calcNodesOther) > 0) && !any(model$getNodeType(calcNodesOther)=="stoch")){
      warning("There are no stochastic nodes in the calcNodesOther provided for Laplace or AGHQuad approximation.")
    }
  }
  if(!calcOtherProvided){
    calcNodesOther <- calcNodesOtherDefault
  }
  if(calcOtherProvided && check) {
    calcOtherCheck <- setdiff(calcNodesOtherDefault, calcNodesOther)
    if(length(calcOtherCheck)) {
      # We only check missing stochastic nodes; determ nodes will be added below
      missingStochNodesInds <- which((model$getNodeType(calcOtherCheck)) == "stoch")
      lengthMissingStochNodes <- length(missingStochNodesInds)
      if(lengthMissingStochNodes){
        missingStochNodes <- calcOtherCheck[missingStochNodesInds]
        errorNodes <- paste0(head(missingStochNodes, n = 4), sep = "", collapse = ", ")
        if(lengthMissingStochNodes > 4) errorNodes <- paste(errorNodes, "...")
        warning(paste0("There are some model nodes (stochastic) that look like they should be\n",
                       "included in the calcNodesOther for parts of the likelihood calculation\n",
                       "outside of Laplace or AGHQuad approximation:\n",
                       errorNodes, "\n",
                       "To silence this warning, include \'check = FALSE\' in the control list\n",
                       "to buildLaplace or as an argument to setupMargNodes."))
      }
    }
    # Check redundant stochastic nodes
    calcOtherCheck <- setdiff(calcNodesOther, calcNodesOtherDefault)
    stochCalcOtherCheck <- calcOtherCheck[model$getNodeType(calcOtherCheck)=="stoch"]
    errorNodes <- stochCalcOtherCheck
    # Check redundant determ nodes
    # N.B. I commented-out this deterministic node checking for reasons similar to above. -Perry
    ## determCalcOtherCheck <- setdiff(calcOtherCheck, stochCalcOtherCheck)
    ## lengthDetermCalcOtherCheck <- length(determCalcOtherCheck)
    ## errorNodes <- character(0)
    ## if(lengthDetermCalcOtherCheck){
    ##   if(!paramDetermDepsCalculated) {
    ##     paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE, includePredictive = FALSE)
    ##     paramDetermDepsCalculated <- TRUE
    ##   }
    ##   for(i in 1:lengthDetermCalcOtherCheck){
    ##     if(!(determCalcOtherCheck[i] %in% paramDetermDeps) ||
    ##        !(any(model$getDependencies(determCalcOtherCheck[i], self = FALSE) %in% calcNodesOtherDefault))){
    ##       errorNodes <- c(errorNodes, determCalcOtherCheck[i])
    ##     }
    ##   }
    ## }
    ## errorNodes <- c(stochCalcOtherCheck, errorNodes)
    if(length(errorNodes)){
      outErrorNodes <- paste0(head(errorNodes, n = 4), sep = "", collapse = ", ")
      if(length(errorNodes) > 4) outErrorNodes <- paste(outErrorNodes, "...")
      warning(paste0("There are some nodes provided in calcNodesOther that look like\n",
                     "they are not needed for parts of the likelihood calculation\n",
                     "outside of Laplace or AGHQuad approximation:\n",
                     outErrorNodes, "\n",
                     "To silence this warning, include \'check = FALSE\' in the control list\n",
                     "to buildLaplace or as an argument to setupMargNodes."))
    }
  }
  # Check and add necessary (upstream) deterministic nodes into calcNodesOther
  # This ensures that deterministic nodes between paramNodes and calcNodesOther are used.
  num_calcNodesOther <- length(calcNodesOther)
  if(num_calcNodesOther > 0){
    if(!paramDetermDepsCalculated) {
      paramDetermDeps <- model$getDependencies(paramNodes, determOnly = TRUE, includePredictive = FALSE)
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
        # givenNodes should only be stochastic
        givenNodes <- setdiff(c(paramNodes, calcNodes),
                              c(randomEffectsNodes,
                                model$getDependencies(randomEffectsNodes, determOnly=TRUE)))
        reSets <- model$getConditionallyIndependentSets(
          nodes = randomEffectsNodes, givenNodes = givenNodes,
          unknownAsGiven = TRUE)
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

## Main function for Adaptive Gauss-Hermite Quadrature
#' @rdname laplace
#' @export
buildAGHQuad <- nimbleFunction(
  name = 'AGHQuad',
  setup = function(model, nQuad = 1, paramNodes, randomEffectsNodes, calcNodes,
                   calcNodesOther, control = list()) {
    split <- extractControlElement(control, 'split', TRUE)
    check <- extractControlElement(control, 'check', TRUE)
    innerOptimWarning <- extractControlElement(control, 'innerOptimWarning', FALSE)

    if(nQuad > 35) {
      print("  Warning: Currently only a maximum of 35 quadrature points are allowed, setting nQuad to 35.")
      nQuad <- 35
    }
    nQuad_ <- nQuad
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
    outerOptimControl_   <- nimOptimDefaultControl()
    innerOptimControl_ <- nimOptimDefaultControl()
    optimControlArgNames <- c("trace", "fnscale", "parscale", "ndeps", "maxit", "abstol", "reltol", "alpha", 
                              "beta", "gamma", "REPORT", "type", "lmm", "factr", "pgtol", "temp", "tmax")
    if(!is.null(control$outerOptimControl)){
      validNames <- intersect(names(control$outerOptimControl), optimControlArgNames)
      numValidNames <- length(validNames)
      if(numValidNames > 0){
        for(i in 1:numValidNames){
          outerOptimControl_[[validNames[i]]] <- control$outerOptimControl[[validNames[i]]]
        }   
      }
    }
    if(!is.null(control$innerOptimControl)) {
      validNames_inner <- intersect(names(control$innerOptimControl), optimControlArgNames)
      numValidNames_inner <- length(validNames_inner)
      if(numValidNames_inner > 0){
        for(i in 1:numValidNames_inner)
          innerOptimControl_[[validNames_inner[i]]] <- control$innerOptimControl[[validNames_inner[i]]]
      }
    }
    outerOptimControl_$fnscale <- -1
    innerOptimControl_$fnscale <- -1
    if(!is.null(control$innerOptimMethod) &&
       ((control$innerOptimMethod %in% c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")) ||
        control$innerOptimMethod %in% ls(nimbleUserNamespace$.optimizers))){
      innerOptimMethod <- control$innerOptimMethod
    }
    else innerOptimMethod <- "BFGS"

    innerOptimStart <- extractControlElement(control, "innerOptimStart", "zero")
    if(!is.character(innerOptimStart) |
       length(innerOptimStart) != 1 |
       !(innerOptimStart %in% (validIOS <- c("last", "last.best", "constant", "random", "model", "zero"))))
      stop(paste("control$innerOptimStart must be one of ", paste0('\'', validIOS, '\'', collapse=",")))

    innerOptimStartValues <- NULL
    if(innerOptimStart == "model") {
      innerOptimStartValues <- 0 # will be ignored but is need to trigger next message
    }
    if(innerOptimStart == "zero") {
      innerOptimStart <- "constant"
      innerOptimStartValues <- 0
    }
    if(!is.null(innerOptimStartValues) &
       !is.null(control$innerOptimStartValues)) {
      message(paste0("ignoring control$innerOptimStartValues because control$innerOptimStart is \"",innerOptimStart,"\""))
    } else {
      innerOptimStartValues <- extractControlElement(control, "innerOptimStartValues", 0)
      if(is.character(innerOptimStartValues))
        if(length(innerOptimStartValues) != 1 |
           !(innerOptimStartValues == "model"))
          stop("The only valid character value for control$innerOptimStartValues is 'model'.")
    }

    ## Create an AGHQuad (Adaptive Gauss-Hermite Quadrature) nimbleFunctionList
    AGHQuad_nfl <- nimbleFunctionList(AGHQuad_BASE)
    scalarRENodes <- model$expandNodeNames(randomEffectsNodes, returnScalarComponents = TRUE)
    nre <- length(scalarRENodes)
    multiSetsCheck <- FALSE ## AGHQ vs Laplace Check in findMLE.
    gridType <- extractControlElement(control, "gridType", "cholesky")
    innerControlList <- list(optimControl=innerOptimControl_,
                             optimMethod=innerOptimMethod,
                             optimStart=innerOptimStart,
                             optimStartValues=innerOptimStartValues,
                             optimWarning=innerOptimWarning,
                             gridType=gridType)
    if(nre > 0){
      ## Record the order of random effects processed internally
      internalRandomEffectsNodes <- NULL
      lenInternalRENodeSets <- NULL
      if(isFALSE(split)) { ## Do all randomEffectsNodes in one set
        internalRandomEffectsNodes <- randomEffectsNodes
        lenInternalRENodeSets <- nre
        reNodesAsScalars <- model$expandNodeNames(internalRandomEffectsNodes, returnScalarComponents = TRUE)

        # old default was "model". new default is "last.best" with values=0
        # Essentially the following steps will now b done in buildOneAGHQuad[1D]
        ## if(is.null(control$innerOptimStart)) innerOptimStart <- values(model, randomEffectsNodes)
        ## else {
        ##   providedStart <- control$innerOptimStart
        ##   if(any(providedStart %in% c("last", "last.best"))) innerOptimStart <- providedStart
        ##   else if(is.numeric(sum(providedStart)) && (length(providedStart) == nre)) innerOptimStart <- providedStart
        ##   else innerOptimStart <- values(model, randomEffectsNodes)
        ## }
        ## ## In case random effects are not properly initialized
        ## if(!any(innerOptimStart %in% c("last", "last.best")) & any(is.infinite(innerOptimStart) | is.na(innerOptimStart) | is.nan(innerOptimStart))){
        ##   all_reTransform <- parameterTransform(model, randomEffectsNodes)
        ##   all_reTransform_length <- all_reTransform$getTransformedLength()
        ##   innerOptimStart <- all_reTransform$inverseTransform(rep(0, all_reTransform_length))
        ## }
        ## Build AGHQuad
        if(nre > 1 | isTRUE(control[['force_nDim']])) {
          AGHQuad_nfl[[1]] <- buildOneAGHQuad(model, nQuad = nQuad_, paramNodes, randomEffectsNodes,
                                              calcNodes,
                                              innerControlList)
          multiSetsCheck <- TRUE
        } else AGHQuad_nfl[[1]] <- buildOneAGHQuad1D(model, nQuad = nQuad_, paramNodes, randomEffectsNodes,
                                                     calcNodes,
                                                     innerControlList)
        }
      else {## Split randomEffectsNodes into conditionally independent sets
        reSets <- MargNodes$randomEffectsSets
        num_reSets <- length(reSets)
        reNodesAsScalars <- character()
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
          these_reNodesAsScalars <- model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE)
          reNodesAsScalars <- c(reNodesAsScalars, these_reNodesAsScalars)
          nre_these <- length(these_reNodesAsScalars)
          lenInternalRENodeSets <- c(lenInternalRENodeSets, nre_these)
          ## Process start values for inner optimisation

          if(is.numeric(innerOptimStartValues) && (length(innerOptimStartValues) == nre)) {
            # input was a vector of all REs, so we must split it accordingly
            if(is.null(names(innerOptimStartValues))) {
              # split by order, because there are no names. use internalRandomEffectsNodes order.
              # the last portion will be for the present set of values
              these_reNodes_inds <- length(reNodesAsScalars) - nre_these + (1:nre_these)
            } else {
              # split by names
              these_reNodes_inds <- match(these_reNodesAsScalars, innerOptimStartValues, nomatch=0)
              if((length(these_reNodes_inds) != nre_these) |
                   (length(unique(these_reNodes_inds)) != length(these_reNodes_inds)) |
                   any(these_reNodes_inds==0))
                warning("There appears to be an incorrect name in control$innerOptimStartValues")
            }
            these_innerOptimStartValues <- innerOptimStartValues[these_reNodes_inds]
          } else
            these_innerOptimStartValues <- innerOptimStartValues

          ## if(is.null(control$innerOptimStart)) innerOptimStart <- values(model, these_reNodes)
          ## else {
          ##   providedStart <- control$innerOptimStart
          ##   if(any(providedStart %in% c("last", "last.best"))) innerOptimStart <- providedStart
          ##   else if(is.numeric(sum(providedStart)) && (length(providedStart) == nre)){
          ##     these_reNodes_inds <- unlist(lapply(model$expandNodeNames(these_reNodes, returnScalarComponents = TRUE), function(x) {which(scalarRENodes == x)}))
          ##     innerOptimStart <- providedStart[these_reNodes_inds]
          ##   }
          ##   else innerOptimStart <- values(model, these_reNodes)
          ## }
          ## ## In case random effects are not properly initialized
          ## if(!any(innerOptimStart %in% c("last", "last.best")) & any(is.infinite(innerOptimStart) | is.na(innerOptimStart) | is.nan(innerOptimStart))){
          ##   these_reTransform <- parameterTransform(model, these_reNodes)
          ##   these_reTransform_length <- these_reTransform$getTransformedLength()
          ##   innerOptimStart <- these_reTransform$inverseTransform(rep(0, these_reTransform_length))
          ## }
          ## Build AGHQuad for each set
          if(nre_these > 1 | isTRUE(control[['force_nDim']])){
            AGHQuad_nfl[[i]] <- buildOneAGHQuad(model, nQuad = nQuad_, paramNodes, these_reNodes, these_calcNodes,
                                               innerControlList)
            #innerOptimControl_, innerOptimMethod, innerOptimStart, these_innerOptimStartValues)
            multiSetsCheck <- TRUE
          }
          else AGHQuad_nfl[[i]] <- buildOneAGHQuad1D(model, nQuad = nQuad_, paramNodes, these_reNodes, these_calcNodes,
                                                     innerControlList)
          #innerOptimControl_, innerOptimMethod, innerOptimStart, these_innerOptimStartValues)
        }
      }
      if(length(lenInternalRENodeSets) == 1) lenInternalRENodeSets <- c(lenInternalRENodeSets, -1)
      reTransform <- parameterTransform(model, internalRandomEffectsNodes)
      reTransform_length <- reTransform$getTransformedLength()
      if(reTransform_length > 1) reTransform_indices <- 1:reTransform_length
      else reTransform_indices <- c(1, -1)
      
      reNodesAsScalars_vec <- reNodesAsScalars
      if(nre == 1) reNodesAsScalars_vec <- c(reNodesAsScalars, "_EXTRA_")
      reNodesAsScalars_first <- reNodesAsScalars[1]
    }
    else{
      ## No random effects
      lenInternalRENodeSets <- numeric(2)
      reTransform <- parameterTransform(model, paramNodes[1], control = list(allowDeterm = FALSE)) ## Won't be needed at all
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
    paramsTransform <- parameterTransform(model, paramNodes, control = list(allowDeterm = FALSE))
    pTransform_length <- paramsTransform$getTransformedLength()
    if(pTransform_length > 1) pTransform_indices <- 1:pTransform_length
    else pTransform_indices <- c(1, -1)
    
    ## Indicator for removing the redundant index -1 in pTransform_indices
    one_time_fixes_done <- FALSE
    ## Default calculation method for AGHQuad
    computeMethod_ <- extractControlElement(control, "computeMethod", 2)

    useInnerCache_ <- extractControlElement(control, "useInnerCache", TRUE)

    ## The nimbleList definitions AGHQuad_params and AGHQuad_summary
    ## have moved to predefined nimbleLists.
  },## End of setup
  run = function(){},
  methods = list(
    getNodeNamesVec = function(returnParams = logical(0, default = TRUE)) {
      one_time_fixes()
      returnType(character(1))
      if(returnParams) return(paramNodesAsScalars_vec)
      else return(reNodesAsScalars_vec)
    },
    getNodeNameSingle = function(returnParams = logical(0, default = TRUE)) {
      returnType(character())
      if(returnParams) return(paramNodesAsScalars_first)
      else return(reNodesAsScalars_first)
    },
    updateSettings = function(innerOptimMethod = character(0, default="NULL"),
                               innerOptimStart = character(0, default="NULL"),
                               innerOptimStartValues = double(1, default=Inf),
                               innerOptimWarning = integer(0, default = -1),
                               useInnerCache = integer(0, default=-1),
                               nQuad = integer(0, default=-1),
                               gridType = character(0, default="NULL"),
                               innerOptimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                               replace_innerOptimControl = logical(0, default=FALSE),
                               outerOptimControl = optimControlNimbleList(default=nimOptimDefaultControl()),
                               replace_outerOptimControl = logical(0, default=FALSE),
                               computeMethod = integer(0, default=-1)
                               ) {
      # checks
      if(innerOptimStart != "NULL") {
        if(innerOptimStart=="zero") {
          stop("innerOptimStart choice 'zero' is not supported in updateSettings. Use innerOptimStart='constant' and innerOptimStartValues=0 to achieve 'zero' behavior.")
        }
        if(innerOptimStart != "last" & innerOptimStart != "last.best" &
             innerOptimStart != "constant" & innerOptimStart != "random" &
             innerOptimStart != "model") {
          stop("Invalid value for innerOptimStart.")
        }
      }
      if(length(innerOptimStartValues) > 1) {
        if(length(innerOptimStartValues) != nre)
          stop("Length of innerOptimStartValues must be 1 or total number of random effects")
      }
      if(nQuad != -1)  {
        if(nQuad < 1) stop("Choose a positive number of grid points.")
        if(nQuad > 35) stop("Currently only a maximum of 35 quadrature points are allowed.")
        threshold <- log(50000) # in text below too
        for(i in seq_along(AGHQuad_nfl)) {
          if(nQuad * log(lenInternalRENodeSets[i]) > threshold) {
            print("Choice of nQuad would yield >50000 nodes for ",lenInternalRENodeSets[i],
                  " integration dimensions in conditionally independent set ", i)
            stop("That is too many nodes.")
          }
        }
      }
      if(computeMethod != -1) {
        if(!any(c(1, 2, 3) == computeMethod))
          stop("computeMethod must be 1, 2, or 3")
      }
      if(gridType != "NULL") {
        if(gridType != "spectral" & gridType != "cholesky")
          stop("gridType must be either cholesky or spectral.")
      }
      # actions
      one_time_fixes()
      if(nQuad != -1) nQuad_ <<- nQuad
      these_initsValues <- innerOptimStartValues
      iStart <- 1
      for(i in seq_along(AGHQuad_nfl)) {
        if(length(innerOptimStartValues) > 1) {
          these_values <- innerOptimStartValues[iStart:(iStart + lenInternalRENodeSets[i] - 1)]
          iStart <- iStart + lenInternalRENodeSets[i]
        }
        AGHQuad_nfl[[i]]$updateSettings(optimMethod = innerOptimMethod,
                                         optimStart = innerOptimStart,
                                         optimStartValues = innerOptimStartValues,
                                         optimWarning = innerOptimWarning,
                                         useInnerCache = useInnerCache,
                                         nQuad = nQuad_,
                                         gridType = gridType,
                                         optimControl = innerOptimControl,
                                         replace_optimControl = replace_innerOptimControl)
      }
      # TO-DO: create useInnerCache_ and allow control arg.
      if(useInnerCache != -1) useInnerCache_ <<- useInnerCache != 0
      if(computeMethod != -1) computeMethod_ <<- computeMethod
      if(replace_outerOptimControl) {
        outerOptimControl$fnscale <- -1
        outerOptimControl_ <<- outerOptimControl
      }
    },
    ## setMethod = function(method = integer()) {
    ##   if(nre == 0) print("AGHQuad or Laplace approximation is not needed for the given model: no random effects")
    ##   if(!any(c(1, 2, 3) == method)) stop("Choose a valid method ID from 1, 2, and 3")
    ##   methodID <<- method
    ## },
    ## getMethod = function() {
    ##   return(methodID)
    ##   returnType(integer())
    ## },
    ## Let the user experiment with different quadrature grids:
    ## setQuadSize = function(nQUpdate = integer()){
    ##   nQuad0 <- nQuad_
    ##   if(nQUpdate < 1) stop("Choose a positive number of grid points.")
    ##   if(nQUpdate > 35) stop("Currently only a maximum of 35 quadrature points are allowed.")
    ##   nQuad_ <<- nQUpdate
    ##   for(i in seq_along(AGHQuad_nfl)) {
    ##     if( lenInternalRENodeSets[i]^nQuad_ > 50000 ){
    ##       nQuad_ <<- nQuad0
    ##       stop("You have exceeded the maximum quadrature grid of 50,000 points.")
    ##     }
    ##     AGHQuad_nfl[[i]]$set_nQuad(nQuad_)
    ##   }
    ## },
    ## setAGHQTransformation = function(method = character()){
    ##   if(method != "spectral" & method != "cholesky") stop("Must choose either cholesky or spectral.")
    ##   for(i in seq_along(AGHQuad_nfl)) AGHQuad_nfl[[i]]$set_transformation(transformation = method)
    ## },
    ## setInnerOptimInits = function(method = character(0), values = double(1)){
    ##   full_values <- length(values) > 1
    ##   if(full_values)
    ##     if(length(values) != nre)
    ##       stop("values may be empty, or have length = 1 or to the total number of scalar random effects.")
    ##   if(length(values) == 0) values <- c(0)
    ##   if(!full_values) these_values <- values
    ##   iStart <- 1
    ##   for(i in seq_along(AGHQuad_nfl)) {
    ##     if(full_values) {
    ##       these_values <- values[ iStart:(iStart + lenInternalRENodeSets[i] - 1) ]
    ##       iStart <- iStart + lenInternalRENodeSets[i]
    ##     }
    ##     AGHQuad_nfl[[i]]$set_reInitMethod(method, these_values)
    ##   }
    ## },
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
    ## Check to see if the inner optimizations converged.
    checkInnerConvergence = function(message = logical(0, default = FALSE)){
      converged <- 0
      for(i in seq_along(AGHQuad_nfl)){
        conCheck <- AGHQuad_nfl[[i]]$check_convergence()
        if(conCheck != 0) {
          converged <- 1
          if(message) print("  [Warning] Inner optimization did not converge for conditionally independent set ", i, " with code ", conCheck, ".")
        }
      }
      returnType(double())
      return(converged)
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
      checkInterrupt()
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
          if(computeMethod_ == 1) ans <- ans + AGHQuad_nfl[[i]]$calcLogLik1(p)
          else if(computeMethod_ == 2) ans <- ans + AGHQuad_nfl[[i]]$calcLogLik2(p)
          else ans <- ans + AGHQuad_nfl[[i]]$calcLogLik3(p)
        }
      }
      if(is.nan(ans) | is.na(ans)) ans <- -Inf
      return(ans)
      returnType(double())
    },
    calcLaplace = function(p = double(1), trans = logical(0, default = FALSE)) {
      if(nQuad_ > 1) {
        stop("Must set nQuad to 1 in order to call calcLaplace. Either call calcLogLik or use updateSettings() to change nQuad.")
      }
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
          if(computeMethod_ == 1) ans <- ans + AGHQuad_nfl[[i]]$gr_logLik1(p)
          else if(computeMethod_ == 2) ans <- ans + AGHQuad_nfl[[i]]$gr_logLik2(p)
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
      if(nQuad_ > 1) {
        stop("Must set nQuad to 1 in order to call calcLaplace. Either call calcLogLik or use updateSettings() to change nQuad.")
      }
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
      cache_outer_logLik(ans) ## Save outer in the inner to cache values at outer mode.
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
    ## Prior contribution to the posterior
		calcPrior_p = function(p = double(1)){
      ## Prior log likelihood:
      values(model, paramNodes) <<- p
      ans <- model$calculate(paramNodes)
      return(ans)
      returnType(double())
		},
    ## Prior contribution to the posterior on the transformed scale.
		calcPrior_pTransformed = function(pTransform = double(1)) {
      p <- paramsTransform$inverseTransform(pTransform)
			ans <- calcPrior_p(p) + logDetJacobian(pTransform)
			return(ans)
			returnType(double())
		},
    ## Calculate posterior density at p log likelihood + log prior.
		calcPostLogDens = function(p = double(1), trans = logical(0, default = FALSE)) {
      ans <- 0
      if(trans) {
        pstar <- paramsTransform$inverseTransform(p)  ## Just want to do this once.
        ans <- ans + logDetJacobian(p)  ## p is transformed, add Jacobian here.
      }else{
        pstar <- p
      }
      ## Error checking when calling calcLogLik.
			ans <- ans + calcLogLik(pstar, FALSE) + calcPrior_p(pstar)
			returnType(double())
			return(ans)
		},
    ## Calculate posterior density at p transformed, log likelihood + log prior (transformed).
		calcPostLogDens_pTransformed = function(pTransform = double(1)) {
      
      ans <- calcPostLogDens(pTransform, TRUE)
      cache_outer_logLik(ans) ## Update internal cache w/ prior.

      if(is.nan(ans) | is.na(ans)) ans <- -Inf			
      returnType(double())
			return(ans)
		},
    ## Gradient of log det jacobian for parameter transformations.
		gr_logDetJacobian = function(pTransform = double(1))
		{
			ans <- derivs(logDetJacobian(pTransform), wrt = pTransform_indices, order = 1)
			return(ans$jacobian[1,])
			returnType(double(1))
		},
    ## Gradient of prior distribution.
		gr_prior = function(p = double(1))
		{
			ans <- derivs(calcPrior_p(p), wrt = p_indices, order = 1)
			return(ans$jacobian[1,])
			returnType(double(1))
		},
    ## Gradient of posterior density on the transformed scale.
		gr_postLogDens_pTransformed = function(pTransform = double(1))
		{
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
    ## For internal purposes of building the gradient
		logDetJacobian = function(pTransform = double(1)){
			ans <- paramsTransform$logDetJacobian(pTransform)
			return(ans)
			returnType(double())
		},
    ## Calculate MLE of parameters    
    findMLE = function(pStart  = double(1, default = Inf),
                       method  = character(0, default = "BFGS"),
                       hessian = logical(0, default = TRUE) ){
      mleRes <- optimize(pStart  = pStart,
                       method  = method,
                       hessian = hessian,
                       parscale = "real") 
      return(mleRes)
      returnType(optimResultNimbleList())
    },
    ## General Maximization Function (Name check: optimize? @perry or Chris?)
    optimize = function(pStart  = double(1, default = Inf),
                       method  = character(0, default = "BFGS"),
                       hessian = logical(0, default = TRUE),
                       parscale = character(0, default = "transformed")) {
      if(!one_time_fixes_done) one_time_fixes() ## Otherwise summary will look bad.
      if(multiSetsCheck & nQuad_ > 1) stop("Currently only Laplace is supported for Maximization. Use setNodeSize(1) to change.")
      if(any(abs(pStart) == Inf)) pStart <- values(model, paramNodes)
      if(length(pStart) != npar) {
        print("  [Warning] For Maximization, pStart should be length ", npar, " but is length ", length(pStart), ".")
        ans <- optimResultNimbleList$new()
        return(ans)
      # stop("Wrong length for pStart in findMLE.")
      }
      ## Reset log likelihood internally for cache.
      reset_outer_inner_logLik()
            
      ## In case parameter nodes are not properly initialized
      if(any_na(pStart) | any_nan(pStart) | any(abs(pStart)==Inf)) pStartTransform <- rep(0, pTransform_length)
      else pStartTransform <- paramsTransform$transform(pStart)
      ## In case bad start values are provided
      if(any_na(pStartTransform) | any_nan(pStartTransform) | any(abs(pStartTransform)==Inf)) pStartTransform <- rep(0, pTransform_length)
     
      optRes <- optim(pStartTransform, calcLogLik_pTransformed, gr_logLik_pTransformed, method = method, control = outerOptimControl_, hessian = hessian)
      
      if(optRes$convergence != 0) 
        print("Warning: optim has a non-zero convergence code: ", optRes$convergence, ".\n",
              "The control parameters of optim can be adjusted in the control argument of\n",
              "buildLaplace or buildAGHQuad via list(outerOptimControl = list()).")
      
      ## Print out warning about inner convergence.
      if( checkInnerConvergence(FALSE) != 0 )
        print("  Warning: inner optimization had a non-zero convergence code. Use checkInnerConvergence(TRUE) to see details.")

      ## Back transform results to original scale if requested.
      p <- paramsTransform$inverseTransform(optRes$par)
      if(parscale == "real") optRes$par <- p
      setModelValues(p) ## Make sure the model object contains all the updated parameter values.
 
      ## Returns on transformed scale just like optim.
      return(optRes)
      returnType(optimResultNimbleList())
    },
    ## User can update whether or not a warning is set for inner optimization.
    ## setInnerOptimWarning = function(warn = logical(0, default = FALSE)){
    ##   for(i in seq_along(AGHQuad_nfl)){
    ##     AGHQuad_nfl[[i]]$set_warning(warn)
    ##   }
    ## },
    ## Grab the inner Cholesky from the cached last values.
    cache_outer_logLik = function(logLikVal = double()){
      for(i in seq_along(AGHQuad_nfl)){
        AGHQuad_nfl[[i]]$save_outer_logLik(logLikVal)
      }
    },
    ## Set cached log lik values to -Inf internally.
    reset_outer_inner_logLik = function(){
      for(i in seq_along(AGHQuad_nfl)){
        AGHQuad_nfl[[i]]$reset_outer_logLik()
      }
    },    
    ## Grab the inner Cholesky from the cached last values.
    get_inner_cholesky = function(atOuterMode = integer(0, default = 0)){
      if(nre == 0) stop("No random effects in the model")
      cholesky <- matrix(value = 0, nrow = nre, ncol = nre)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        numre <- lenInternalRENodeSets[i]
        cholesky[(tot+1):(tot+numre), (tot+1):(tot+numre)] <- AGHQuad_nfl[[i]]$get_inner_negHessian_chol(atOuterMode)
        tot <- tot + numre
      }
      return(cholesky)
      returnType(double(2))
    },
    ## Grab the inner mode from the cached last values.
    get_inner_mode = function(atOuterMode = integer(0, default = 0)){
      if(nre == 0) stop("No random effects in the model")
      raneff <- numeric(nre)
      tot <- 0
      for(i in seq_along(AGHQuad_nfl)){
        numre <- lenInternalRENodeSets[i]
        raneff[(tot+1):(tot+numre)] <- AGHQuad_nfl[[i]]$get_inner_mode(atOuterMode)
        tot <- tot + numre
      }
      return(raneff)
      returnType(double(1))
    },    
    ## Optimized random effects given transformed parameter values
    optimRandomEffects = function(pTransform = double(1)){
      if(nre == 0) stop("No random effects in the model")
      p <- pInverseTransform(pTransform)
      raneff <- numeric(nre)
      tmp <- numeric(nre) ## Not sure this is needed. 
      tot <- 0

      computeMethod <- -1
      if(useInnerCache_){
        pMLE <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 1)
        pLast <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 0)
        ## Cache check for either last value or MLE
        if(all(p == pMLE)) computeMethod <- 1
        else if(all(p == pLast)) computeMethod <- 0
      }

      for(i in seq_along(AGHQuad_nfl)){
        if(computeMethod == -1 ){
          if(computeMethod_ == 1) tmp <- AGHQuad_nfl[[i]]$update_max_inner_logLik_internal(p)
          else tmp <- AGHQuad_nfl[[i]]$update_max_inner_logLik(p)
        }else{
          tmp <- AGHQuad_nfl[[i]]$get_inner_mode(atOuterMode = computeMethod)
        }
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

      outer_mode_case <- -1
      if(useInnerCache_){
        pMLE <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 1)
        pLast <- AGHQuad_nfl[[1]]$get_param_value(atOuterMode = 0)
        ## Cache check for either last value or MLE
        if(all(p == pMLE)) outer_mode_case <- 1
        else if(all(p == pLast)) outer_mode_case <- 0
      }

      for(i in seq_along(AGHQuad_nfl)){
        numre <- lenInternalRENodeSets[i]
        if(outer_mode_case == -1){
          tmp <- AGHQuad_nfl[[i]]$negHess(p, reTransform[(tot+1):(tot+numre)])
        }else{
          U <- AGHQuad_nfl[[i]]$get_inner_negHessian_chol(atOuterMode = outer_mode_case)
          tmp <- t(U) %*% U
        }
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
        if(computeMethod_ == 1) tmp <- AGHQuad_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform[(tot+1):(tot+numre)])
        else tmp <- AGHQuad_nfl[[i]]$hess_joint_logLik_wrt_p_wrt_re(p, reTransform[(tot+1):(tot+numre)])
        ans[1:npar, (tot+1):(tot+numre)] <- tmp
        tot <- tot + numre
      }
      return(ans)
      returnType(double(2))
    },
    ## Gives the user control to start fresh by removing internally saved values.
    ## setInnerCache = function(useCache = logical(0, default = TRUE)){
    ##   innerCache <<- useCache
    ##   for(i in seq_along(AGHQuad_nfl)) AGHQuad_nfl[[i]]$set_inner_cache(useCache)
    ## },
    ## Set all model values after finding the MLE. Function will repeat inner optimization if the inner cached values
    ## the inner cached values don't match p.
    setModelValues = function(p = double(1)){
      for(i in seq_along(AGHQuad_nfl))
         AGHQuad_nfl[[i]]$set_randomeffect_values(p)
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
        optreTransform <- optimRandomEffects(pTransform)  ## *** Replace this with cached inner modes.
        optre <- reInverseTransform(optreTransform)
        ntot <- npar + nre
        if(jointCovariance) {
          ## Inverse of the negative Hessian of log-likelihood wrt transformed random effects at MLEs
          inv_negHess <- inverse_negHess(p, optreTransform)   ## *** Replace this with cached inner modes.
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
										 
#' Summarize results from Laplace or adaptive Gauss-Hermite quadrature approximation
#'
#' Process the results of the `findMLE` method of a nimble Laplace or AGHQ approximation
#' into a more useful format.
#'
#' @param laplace The Laplace approximation object, typically the compiled one.
#'   This would be the result of compiling an object returned from
#'   `buildLaplace`.
#'
#' @param AGHQ Same as \code{laplace}. Note that `buildLaplace` and
#'   `buildAGHQuad` create the same kind of algorithm object that can be used
#'   interchangeably. `buildLaplace` simply sets the number of quadrature points
#'   (`nQuad`) to 1 to achieve Laplace approximation as a special case of AGHQ.
#'
#' @param MLEoutput The maximum likelihood estimate using Laplace or AGHQ,
#'   returned from e.g. `approx$findMLE(...)`, where \code{approx} is the
#'   algorithm object returned by `buildLaplace` or `buildAGHQuad`, or (more
#'   typically) the result of compiling that object with `compileNimble`. See
#'   `help(buildLaplace)` for more information.
#'
#' @param originalScale Should results be returned using the original
#'   parameterization in the model code (TRUE) or the potentially transformed
#'   parameterization used internally by the Laplace approximation (FALSE).
#'   Transformations are used for any parameters and/or random effects that have
#'   constrained ranges of valid values, so that in the transformed parameter
#'   space there are no constraints. 
#'
#' @param randomEffectsStdError If TRUE, calculate the standard error of the
#'   estimates of random effects values.
#'
#' @param jointCovariance If TRUE, calculate the joint covariance matrix of
#'   the parameters and random effects together. If FALSE, calculate the 
#'   covariance matrix of the parameters.
#'
#' @details
#'
#' The numbers obtained by this function can be obtained more directly by
#' `approx$summary(...)`. The added benefit of `summaryLaplace` is to arrange
#' the results into data frames (for parameters and random effects), with row
#' names for the model nodes, and also adding row and column names to the
#' covariance matrix.
#'
#' @return
#'
#' A list with data frames `params` and `randomEffects`, each with columns for
#' `estimate` and (possibly) `se` (standard error) and row names for model
#' nodes, a matrix `vcov` with the covariance matrix with row and column names,
#' and `originalScale` with the input value of `originalScale` so it is recorded
#' for later use if wanted.
#'
#' @aliases summaryAGHQ
#'
#' @name summaryLaplace
#'
#' @export
summaryLaplace <- function(laplace, MLEoutput,
                           originalScale =TRUE,
                           randomEffectsStdError = FALSE,
                           jointCovariance = FALSE) {
  summary <- laplace$summary(MLEoutput, originalScale = originalScale,
                             randomEffectsStdError = randomEffectsStdError,
                             jointCovariance = jointCovariance)
  paramNames <- summary$params$names
  paramEsts <- summary$params$estimates
  if(length(paramEsts) < length(paramNames)) paramNames <- paramNames[1:(length(paramNames)-1)]
  names(paramEsts) <- paramNames
  stdErrParams <- summary$params$stdErrors
  paramsDF <- data.frame(estimate = paramEsts, se = stdErrParams, row.names = paramNames)

  REnames <- summary$randomEffects$names
  REests <- summary$randomEffects$estimates
  if(length(REests) < length(REnames)) REnames <- REnames[1:(length(REnames)-1)]
  REstdErrs <- summary$randomEffects$stdErrors
  if(length(REstdErrs))
    REDF <- data.frame(estimate = REests, se = REstdErrs, row.names = REnames)
  else
    REDF <- data.frame(estimate = REests, row.names = REnames)

  vcov <- summary$vcov
  if (dim(vcov)[1] == length(paramNames)) {
      colnames(vcov) <- rownames(vcov) <- c(paramNames)
  } else {
      colnames(vcov) <- rownames(vcov) <- c(paramNames, REnames)
  }
  list(params = paramsDF,
       randomEffects = REDF,
       vcov = vcov,
       originalScale = originalScale)
}

#' @rdname summaryLaplace
#' @export
summaryAGHQ <- function(AGHQ, MLEoutput,
                        originalScale =TRUE,
                        randomEffectsStdError = FALSE,
                        jointCovariance = FALSE) {
  summaryLaplace(AGHQ, MLEoutput, originalScale, randomEffectsStdError, jointCovariance)
}

#' Combine steps of running Laplace or adaptive Gauss-Hermite quadrature approximation
#'
#' Use an approximation (compiled or uncompiled) returned from
#' `buildLaplace` or `buildAGHQuad` to find the maximum likelihood estimate and return it
#' with random effects estimates and/or standard errors.
#'
#' @aliases runAGHQ runLaplace
#'
#' @param laplace A (compiled or uncompiled) nimble laplace approximation object
#'   returned from `buildLaplace` or `buildAGHQuad`. These return the same type of
#'   approximation algorithm object. `buildLaplace` is simply `buildAGHQuad`
#'   with `nQuad=1`.
#'
#' @param AGHQ Same as \code{laplace}.
#'
#' @method Optimization method for outer optimization. See \code{method}
#'   argument to \code{findMLE} method in \link{\code{buildLaplace}}.
#'
#' @param pStart Initial values for parameters to begin optimization search for
#'   the maximum likelihood estimates. If omitted, the values currently in the
#'   (compiled or uncompiled) model object will be used.
#'
#' @param originalScale If \code{TRUE}, return all results on the original scale
#'   of the parameters and/or random effects as written in the model. Otherwise,
#'   return all results on potentially unconstrained transformed scales that are
#'   used in the actual computations. Transformed scales (parameterizations) are
#'   used if any parameter or random effect has contraint(s) on its support
#'   (range of allowed values). Default = \code{TRUE}.
#'
#' @param randomEffectsStdError If \code{TRUE}, include standard errors for the
#'   random effects estimates. Default = \code{TRUE}.
#'
#' @param jointCovariance If \code{TRUE}, return the full joint covariance
#'   matrix (inverse of the Hessian) of parameters and random effects. Default =
#'   \code{FALSE}.
#'
#' @details
#'
#' Adaptive Gauss-Hermite quadrature is a generalization of Laplace
#' approximation. \code{runLaplace} simply calles \code{runAGHQ} and provides a
#' convenient name.
#'
#' These functions manage the steps of calling the `findMLE` method to obtain
#' the maximum likelihood estimate of the parameters and then the
#' `summaryLaplace` function to obtain standard errors, (optionally) random
#' effects estimates (conditional modes), their standard errors, and the full
#' parameter-random effects covariance matrix.
#'
#' Note that for `nQuad > 1` (see \link{\code{buildAGHQuad}}), i.e., AGHQ with
#' higher order than Laplace approximation, maximum likelihood estimation is
#' available only if all random effects integrations are univariate. With
#' multivariate random effects integrations, one can use `nQuad > 1` only to
#' calculate marginal log likelihoods at given parameter values. This is useful
#' for checking the accuracy of the log likelihood at the MLE obtained for
#' Laplace approximation (`nQuad == 1`). `nQuad` can be changed using the
#' `updateSettings` method of the approximation object.
#'
#' See \link{\code{summaryLaplace}}, which is called for the summary components.
#'
#' @return
#'
#' A list with elements \code{MLE} and \code{summary}.
#'
#' \code{MLE} is the result of the \code{findMLE} method, which contains the
#' parameter estimates and Hessian matrix. This is considered raw output, and
#' one should normally use instead the contents of \code{summary}. (For example
#' not that the Hessian matrix in \code{MLE} may not correspond to the same
#' scale as the parameter estimates if a transformation was used to operate in
#' an unconstrained parameter space.)
#'
#' \code{summary} is the result of \code{summaryLaplace} (or equivalently
#' \code{summaryAGHQ}), which contains parameter estimates and standard errors,
#' and optionally other requested components. All results in this object will be
#' on the same scale (parameterization), either original or transformed, as
#' requested.
#'
#' @export
runLaplace <- function(laplace, ...) {
  runAGHQ(AGHQ = laplace, ...)
}

#' @rdname runLaplace
#' @export
runAGHQ <- function(AGHQ, pStart, method = "BFGS",
                    originalScale = TRUE,
                    randomEffectsStdError = TRUE,
                    jointCovariance = FALSE) {
  if(missing(AGHQ)) stop('must provide a NIMBLE Laplace or AGHQ algorithm')
  if(!identical(nfGetDefVar(AGHQ, 'name'), 'AGHQuad'))
    stop('AGHQ or laplace argument must be a NIMBLE Lapace orAGHQ algorithm (compiled or uncompiled) from buildLaplace or buildAGHQuad.')
  if(!is.Cnf(AGHQ))
    messageIfVerbose('  [Warning] Running an uncompiled Laplace or AGHQ algorithm.',
                     ' Use compileNimble() for faster execution.')

  if(missing(pStart)) pStart <- Inf # code to use values in the model

  opt <- try(AGHQ$findMLE(pStart = pStart, method = method, hessian = TRUE))
  if(inherits(opt, "try-error"))
    stop("method findMLE had an error.")

  summary  <- try(summaryLaplace(laplace=AGHQ, MLEoutput=opt,
                              originalScale=originalScale,
                              randomEffectsStdError=randomEffectsStdError,
                              jointCovariance=jointCovariance))
  if(inherits(summary, "try-error")) {
    warning("summaryLaplace had an error. Only the MLE result will be returned.")
    summary <- NULL
  }
  list(MLE = opt, summary=summary)
}

#' Laplace approximation and adaptive Gauss-Hermite quadrature
#' 
#' Build a Laplace or AGHQ approximation algorithm for a given NIMBLE model.
#' 
#' @param model a NIMBLE model object, such as returned by \code{nimbleModel}.
#'   The model must have automatic derivatives (AD) turned on, e.g. by using
#'   \code{buildDerivs=TRUE} in \code{nimbleModel}.
#' @param nQuad number of quadrature points for AGHQ (in one dimension). Laplace approximation is
#'   AGHQ with `nQuad=1`. Only odd numbers of nodes really
#'   make sense. Often only one or a few nodes can achieve high accuracy. A maximum of
#'   35 nodes is supported. Note that for multivariate quadratures, the number
#'   of nodes will be (number of dimensions)^nQuad.
#' @param paramNodes a character vector of names of parameter nodes in the
#'   model; defaults are provided by \code{\link{setupMargNodes}}.
#'   Alternatively, \code{paramNodes} can be a list in the format returned by
#'   \code{setupMargNodes}, in which case \code{randomEffectsNodes},
#'   \code{calcNodes}, and \code{calcNodesOther} are not needed (and will be
#'   ignored).
#' @param randomEffectsNodes a character vector of names of continuous
#'   unobserved (latent) nodes to marginalize (integrate) over using Laplace
#'   approximation; defaults are provided by \code{\link{setupMargNodes}}.
#' @param calcNodes a character vector of names of nodes for calculating the
#'   integrand for Laplace approximation; defaults are provided by
#'   \code{\link{setupMargNodes}}. There may be deterministic nodes between
#'   \code{paramNodes} and \code{calcNodes}. These will be included in
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
#' @param control a named list for providing additional settings used in Laplace
#'   approximation. See \code{control} section below. Most of these can be
#'   updated later with the `updateSettings` method.
#'
#' @section \code{buildLaplace}:
#'
#' \code{buildLaplace} creates an object that can run Laplace approximation and
#'   for a given model or part of a model. \code{buildAGHQuad} creates an object
#'   that can run adaptive Gauss-Hermite quadrature (AGHQ, sometimes called
#'   "adaptive Gaussian quadrature") for a given model or part of a model.
#'   Laplace approximation is AGHQ with one quadrature point, hence
#'   `buildLaplace` simply calls `buildAGHQ` with `nQuad=1`. These methods
#'   approximate the integration over continuous random effects in a
#'   hierarchical model to calculate the (marginal) likelihood.
#'
#' \code{buildAGHQ} and \code{buildLaplace} will by default (unless changed
#' manually via `control$split`) determine from the model which random effects
#' can be integrated over (marginalized) independently. For example, in a GLMM
#' with a grouping factor and an independent random effect intercept for each
#' group, the random effects can be marginalized as a set of univariate
#' approximations rather than one multivariate approximation. On the other hand,
#' correlated or nested random effects would require multivariate marginalization.
#'
#' Maximum likelihood estimation is available for Laplace approximation
#' (`nQuad=1`) with univariate or multivariate integrations. With `nQuad > 1`,
#' maximum likelihood estimation is available only if all integrations are
#' univariate (e.g., a set of univariate random effects). If there are
#' multivariate integrations, these can be calculated at chosen input parameters
#' but not maximized over parameters. For example, one can find the MLE based on
#' Laplace approximation and then increase `nQuad` (using the `updateSettings`
#' method below) to check on accuracy of the marginal log likelihood at the MLE.
#'
#' Beware that quadrature will use `nQuad^k` quadrature points, where `k` is the
#' dimension of each integration. Therefore quadrature for `k` greater that 2 or
#' 3 can be slow. As just noted, `buildAGHQuad` will determine independent
#' dimensions of quadrature, so it is fine to have a set of univariate random
#' effects, as these will each have k=1. Multivariate quadrature (k>1) is only
#' necessary for nested, correlated, or otherwise dependent random effects.
#'
#' The recommended way to find the maximum likelihood estimate and associated
#' outputs is by calling \link{\code{runLaplace}} or \link{\code{runAGHQ}}. The
#' input should be the compiled Laplace or AGHQ algorithm object. This would be
#' produced by running \link{\code{compileNimble}} with input that is the result
#' of \code{buildLaplace} or \code{buildAGHQuad}.
#'
#' For more granular control, see below for methods \code{findMLE} and
#'   \code{summary}. See function \link{\code{summaryLaplace}} for an easier way
#'   to call the \code{summary} method and obtain results that include node
#'   names. These steps are all done within \code{runLaplace} and
#'   \code{runAGHQ}.
#'
#' The NIMBLE User Manual at r-nimble.org also contains an example of Laplace
#' approximation.
#'
#' @section How input nodes are processed
#'
#' \code{buildLaplace} and \code{buildAGHQuad} make good tries at deciding what
#' to do with the input model and any (optional) of the node arguments. However,
#' random effects (over which approximate integration will be done) can be
#' written in models in multiple equivalent ways, and customized use cases may
#' call for integrating over chosen parts of a model. Hence, one can take full
#' charge of how different parts of the model will be used.
#'
#' Any of the input node vectors, when provided, will be processed using
#'   \code{nodes <- model$expandNodeNames(nodes)}, where \code{nodes} may be
#'   \code{paramNodes}, \code{randomEffectsNodes}, and so on. This step allows
#'   any of the inputs to include node-name-like syntax that might contain
#'   multiple nodes. For example, \code{paramNodes = 'beta[1:10]'} can be
#'   provided if there are actually 10 scalar parameters, 'beta[1]' through
#'   'beta[10]'. The actual node names in the model will be determined by the
#'   \code{exapndNodeNames} step.
#'
#' In many (but not all) cases, one only needs to provide a NIMBLE model object
#'   and then the function will construct reasonable defaults necessary for
#'   Laplace approximation to marginalize over all continuous latent states 
#'   (aka random effects) in a model. The default values for the four groups of 
#'   nodes are obtained by calling \code{\link{setupMargNodes}}, whose arguments 
#'   match those here (except for a few arguments which are taken from control 
#'   list elements here).
#'
#' \code{setupMargNodes} tries to give sensible defaults from
#'   any combination of \code{paramNodes}, \code{randomEffectsNodes},
#'   \code{calcNodes}, and \code{calcNodesOther} that are provided. For example,
#'   if you provide only \code{randomEffectsNodes} (perhaps you want to
#'   marginalize over only some of the random effects in your model),
#'   \code{setupMargNodes} will try to determine appropriate choices for the
#'   others.
#'
#' \code{setupMargNodes} also determines which integration dimensions are
#' conditionally independent, i.e., which can be done separately from each
#' other. For example, when possible, 10 univariate random effects will be split
#' into 10 univariate integration problems rather than one 10-dimensional
#' integration problem.
#'
#' The defaults make general assumptions such as that
#'   \code{randomEffectsNodes} have \code{paramNodes} as parents. However, The
#'   steps for determining defaults are not simple, and it is possible that they
#'   will be refined in the future. It is also possible that they simply don't
#'   give what you want for a particular model. One example where they will not
#'   give desired results can occur when random effects have no prior
#'   parameters, such as `N(0,1)` nodes that will be multiplied by a scale
#'   factor (e.g. sigma) and added to other explanatory terms in a model. Such
#'   nodes look like top-level parameters in terms of model structure, so
#'   you must provide a \code{randomEffectsNodes} argument to indicate which
#'   they are.
#'
#' It can be helpful to call \code{setupMargNodes} directly to see exactly how
#'   nodes will be arranged for Laplace approximation. For example, you may want
#'   to verify the choice of \code{randomEffectsNodes} or get the order of
#'   parameters it has established to use for making sense of the MLE and
#'   results from the \code{summary} method. One can also call
#'   \code{setupMargNodes}, customize the returned list, and then provide that
#'   to \code{buildLaplace} as \code{paramNodes}. In that case,
#'   \code{setupMargNodes} will not be called (again) by \code{buildLaplace}.
#'
#' If \code{setupMargNodes} is emitting an unnecessary warning, simply use
#'   \code{control=list(check=FALSE)}.
#'
#' @section Managing parameter transformations that may be used internally
#'
#' If any \code{paramNodes} (parameters) or \code{randomEffectsNodes} (random
#'   effects / latent states) have constraints on the range of valid values
#'   (because of the distribution they follow), they will be used on a
#'   transformed scale determined by \code{parameterTransform}. This means the
#'   Laplace approximation itself will be done on the transformed scale for
#'   random effects and finding the MLE will be done on the transformed scale
#'   for parameters. For parameters, prior distributions are not included in
#'   calculations, but they are used to determine valid parameter ranges and
#'   hence to set up any transformations. For example, if \code{sigma} is a
#'   standard deviation, you can declare it with a prior such as \code{sigma ~
#'   dhalfflat()} to indicate that it must be greater than 0.
#'
#' For default determination of when transformations are needed, all parameters
#'   must have a prior distribution simply to indicate the range of valid
#'   values. For a param \code{p} that has no constraint, a simple choice is
#'   \code{p ~ dflat()}.
#'
#' @section Understanding inner and outer optimizations
#'
#' Note that there are two numerical optimizations when finding maximum
#' likelihood estimates with a Laplace or (1D) AGHQ algorithm: (1) maximizing
#' the joint log-likelihood of random effects and data given a parameter value
#' to construct the approximation to the marginal log-likelihood at the given
#' parameter value; (2) maximizing the approximation to the marginal
#' log-likelihood over the parameters. In what follows, the prefix 'inner'
#' refers to optimization (1) and 'outer' refers to optimization (2). Currently
#' both optimizations default to using method \code{"BFGS"}. However, one can
#' use other optimizers or simply run optimization (2) manually from R; see the
#' example below. In some problems, choice of inner and/or outer optimizer can
#' make a big difference for obtaining accurate results, especially for standard
#' errors. Hence it is worth experimenting if one is in doubt.
#'
#' @section \code{control} list arguments
#'
#' The \code{control} list allows additional settings to be made using named
#' elements of the list. Most (or all) of these can be updated later using the
#' `updateSettings` method. Supported elements include:
#'
#' \itemize{
#'
#'   \item \code{split}. If TRUE (default), \code{randomEffectsNodes} will be
#'         split into conditionally independent sets if possible. This
#'         facilitates more efficient Laplace or AGHQ approximation because each
#'         conditionally independent set can be marginalized independently. If
#'         FALSE, \code{randomEffectsNodes} will be handled as one multivariate
#'         block, with one multivariate approximation. If \code{split} is a
#'         numeric vector, \code{randomEffectsNodes} will be split by calling
#'         \code{split}(\code{randomEffectsNodes}, \code{control$split}). The
#'         last option allows arbitrary control over how
#'         \code{randomEffectsNodes} are blocked.
#'
#'   \item \code{check}. If TRUE (default), a warning is issued if
#'         \code{paramNodes}, \code{randomEffectsNodes} and/or \code{calcNodes}
#'         are provided but seem to have missing or unnecessary
#'         elements based on some default inspections of the model. If
#'         unnecessary warnings are emitted, simply set \code{check=FALSE}.
#'
#'   \item \code{innerOptimControl}. A list (either an R list or a
#'         `optimControlNimbleList`) of control parameters for the inner
#'         optimization of Laplace approximation using \code{nimOptim}. See
#'         'Details' of \code{\link{nimOptim}} for further information. Default
#'         is `nimOptimDefaultControl()`.
#'
#'   \item \code{innerOptimMethod}. Optimization method to be used in
#'         \code{nimOptim} for the inner optimization. See 'Details' of
#'         \code{\link{nimOptim}}. Currently \code{nimOptim} in NIMBLE supports:
#'         "\code{Nelder-Mead}", "\code{BFGS}", "\code{CG}", "\code{L-BFGS-B}",
#'         "\code{nlminb}", and user-provided optimizers. By default, method
#'         "\code{BFGS}" is used for both univariate and multivariate cases. For
#'         \code{"nlminb"} or user-provided optimizers, only a subset of
#'         elements of the \code{innerOptimControlList} are supported. (Note
#'         that control over the outer optimization method is available as an
#'         argument to `findMLE`). Choice of optimizers can be important and so
#'         can be worth exploring.
#'
#'   \item \code{innerOptimStart}. Method for determining starting values for
#'         the inner optimization. Options are:
#'
#' \itemize{
#'
#' \item \code{"zero"} (default): use all zeros;
#'
#' \item \code{"last"}: use the result of the last inner optimization;
#'
#' \item \code{"last.best"}: use the result of the best inner
#'         optimization so far for each conditionally independent part of the
#'         approximation;
#'
#' \item \code{"constant"}: always use the same values, determined by
#'         \code{innerOptimStartValues};
#'
#' \item \code{"random"}: randomly draw new starting values from the
#'       model (i.e., from the prior);
#'
#' \item \code{"model"}: use values for random effects stored in the
#'          model, which are determined from the first call.
#'
#' }
#'
#'       Note that \code{"model"} and \code{"zero"} are shorthand for
#'         \code{"constant"} with particular choices of
#'         \code{innerOptimStartValues}. Note that \code{"last"} and
#'         \code{"last.best"} require a choice for the very first values, which will
#'         come from \code{innerOptimStartValues}. The default is
#'         \code{innerOptimStart="zero"} and may change in the future.
#'
#'   \item \code{innerOptimStartValues}. Values for some of
#'         \code{innerOptimStart} approaches. If a scalar is provided, that
#'         value is used for all elements of random effects for each
#'         conditionally independent set. If a vector is provided, it must be
#'         the length of *all* random effects. If these are named (by node
#'         names), the names will be used to split them correctly among each
#'         conditionally independent set of random effects. If they are not
#'         named, it is not always obvious what the order should be because it
#'         may depend on the conditionally independent sets of random
#'         effects. It should match the order of names returned as part of
#'         `summaryLaplace`.
#'
#'   \item \code{innerOptimWarning}. If FALSE (default), do not emit warnings
#'   from the inner optimization. Optimization methods may sometimes emit a
#'   warning such as for bad parameter values encountered during the
#'   optimization search. Often, a method can recover and still find the
#'   optimum. In the approximations here, sometimes the inner optimization
#'   search can fail entirely, yet the outer optimization see this as one failed
#'   parameter value and can recover. Hence, it is often desirable to silence
#'   warnings from the inner optimizer, and this is done by default. Set
#'   \code{innerOptimWarning=TRUE} to see all warnings.
#'
#'   \item \code{useInnerCache}. If TRUE (default), use caching system for
#'     efficiency of inner optimizations. The caching system records one set of
#'     previous parameters and uses the corresponding results if those parameters
#'    are used again (e.g., in a gradient call). This should generally not be
#'    modified.
#'
#'   \item \code{outerOptimControl}. A list of control parameters for maximizing
#'         the Laplace log-likelihood using \code{nimOptim}. See 'Details' of
#'         \code{\link{nimOptim}} for further information.
#'
#'   \item \code{computeMethod}. There are three approaches available for
#'   internal details of how the approximations, and specifically derivatives
#'   involved in their calculation, are handled. These are labeled simply 1, 2,
#'   and 3, and the default is 2. The relative performance of the methods will
#'   depend on the specific model. Users wanting to explore efficiency can try
#'   switching from method 2 (default) to methods 1 or 3 and comparing
#'   performance. The first Laplace approximation with each method will be
#'   (much) slower than subsequent Laplace approximations. Further details are
#'   not provided at this time.
#'
#'  \item \code{gridType} (relevant only \code{nQuad>1}). For multivariate AGHQ,
#'  a grid must be constructed based on the Hessian at the inner mode. Options
#'  include "cholesky" (default) and "spectral" (i.e., eigenvectors and
#'  eigenvalues) for the corresponding matrix decompositions on which the grid
#'  can be based.
#'
#' } # end itemize
#'
#' @section Available methods
#'
#' The object returned by \code{buildLaplace} is a nimbleFunction object with
#' numerous methods (functions). Here these are described in three tiers of user
#' relevance.
#'
#' @section Most useful methods
#'
#' The most relevant methods to a user are:
#'
#' \itemize{
#'
#' \item \code{calcLogLik(p, trans=FALSE)}. Calculate the approximation to the
#'       marginal log-likelihood function at parameter value \code{p}, which (if
#'       \code{trans} is FALSE) should match the order of \code{paramNodes}. For
#'       any non-scalar nodes in \code{paramNodes}, the order within the node is
#'       column-major. The order of names can be obtained from method
#'       \code{getNodeNamesVec(TRUE)}. Return value is the scalar (approximate,
#'       marginal) log likelihood.
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
#' \item \code{calcLaplace(p, trans)}. This is the same as \code{calcLogLik} but
#'        requires that the approximation be Laplace (i.e \code{nQuad} is 1),
#'        and results in an error otherwise.
#'
#' \item \code{findMLE(pStart, method, hessian)}. Find the maximum likelihood
#'         estimates of parameters using the approximated marginal likelihood.
#'         This can be used if \code{nQuad} is 1 (Laplace case) or if
#'         \code{nQuad>1} and all marginalizations involve only univariate
#'         random effects. Arguments include \code{pStart}: initial parameter
#'         values (defaults to parameter values currently in the model);
#'         \code{method}: (outer) optimization method to use in \code{nimOptim}
#'         (defaults to "BFGS", although some problems may benefit from other
#'         choices); and \code{hessian}: whether to calculate and return the
#'         Hessian matrix (defaults to \code{TRUE}, which is required for
#'         subsequent use of `summary` method). Second derivatives in the
#'         Hessian are determined by finite differences of the gradients
#'         obtained by automatic differentiation (AD). Return value is a
#'         nimbleList of type \code{optimResultNimbleList}, similar to what is
#'         returned by R's optim. See \code{help(nimOptim)}. Note that
#'         parameters (`par`) are returned for the natural parameters, i.e. how
#'         they are defined in the model. But the `hessian`, if requested, is
#'         computed for the parameters as transformed for optimization if
#'         necessary. Hence one must be careful interpreting `hessian` if any
#'         parameters have constraints, and the safest next step is to use the
#'         `summary` method or `summaryLaplace` function.
#'
#' \item \code{summary(MLEoutput, originalScale, randomEffectsStdError,
#'        jointCovariance)}. Summarize the maximum likelihood estimation
#'        results, given object \code{MLEoutput} that was returned by
#'        \code{findMLE}. The summary can include a covariance matrix for the
#'        parameters, the random effects, or both), and these can be returned on
#'        the original parameter scale or on the (potentially) transformed
#'        scale(s) used in estimation. It is often preferred instead to call
#'        function (not method) `summaryLaplace` because this will attach
#'        parameter and random effects names (i.e., node names) to the results.
#'
#' In more detail, \code{summary} accepts the following optional arguments:
#'
#'        \itemize{
#'
#'           \item \code{originalScale}. Logical. If TRUE, the function returns
#'           results on the original scale(s) of parameters and random effects;
#'           otherwise, it returns results on the transformed scale(s). If there
#'           are no constraints, the two scales are identical. Defaults to TRUE.
#'
#'           \item \code{randomEffectsStdError}. Logical. If TRUE, standard
#'           errors of random effects will be calculated.
#'           Defaults to FALSE.
#'
#'           \item \code{jointCovariance}. Logical. If TRUE, the joint
#'           variance-covariance matrix of the parameters and the random effects
#'           will be returned. If FALSE, the variance-covariance matrix of the 
#'           parameters will be returned. Defaults to FALSE.
#'
#'        }
#'
#'        The object returned by \code{summary} is an \code{AGHQuad_summary}
#'        nimbleList with elements:
#'
#'        \itemize{
#'
#'           \item \code{params}. A nimbleList that contains estimates and
#'           standard errors of parameters (on the original or transformed
#'           scale, as chosen by \code{originalScale}).
#'
#'           \item \code{randomEffects}. A nimbleList that contains estimates of
#'           random effects and, if requested
#'           (\code{randomEffectsStdError=TRUE}) their standard errors, on
#'           original or transformed scale. Standard errors are calculated
#'           following the generalized delta method of Kass and Steffey (1989).
#'
#'           \item \code{vcov}. If requested (i.e.
#'           \code{jointCovariance=TRUE}), the joint variance-covariance
#'           matrix of the parameters and random effects, on original or
#'           transformed scale. If \code{jointCovariance=FALSE}, the
#'           covariance matrix of the parameters, on original or transformed 
#'           scale.
#'
#'           \item \code{scale}. \code{"original"} or \code{"transformed"}, the
#'           scale on which results were requested.
#'        }
#'     }
#'
#'
#' @section Methods for more advanced uses
#'
#' Additional methods to access or control more details of the Laplace
#' approximation include:
#'
#' \itemize{
#'
#'   \item \code{updateSettings}. This provides a single function through which
#'   many of the settings described above (mostly for the \code{control} list)
#'   can be later changed. Options that can be changed include:
#'   \code{innerOptimMethod}, \code{innerOptimStart},
#'   \code{innerOptimStartValues}, \code{useInnerCache}, \code{nQuad},
#'   \code{gridType}, \code{innerOptimControl}, \code{outerOptimControl}, and
#'   \code{computeMethod}. For \code{innerOptimStart}, method "zero" cannot be
#'   specified but can be achieved by choosing method "constant" with
#'   \code{innerOptimStartValues=0}. Only provided options will be modified. The
#'   exceptions are \code{innerOptimControl}, \code{outerOptimControl}, which
#'   are replaced only \code{replace_innerOptimControl=TRUE} or
#'   \code{replace_outerOptimControl=TRUE}, respectively.
#'
#'   \item \code{getNodeNamesVec(returnParams)}. Return a vector (>1) of names
#'   of parameters/random effects nodes, according to \code{returnParams =
#'   TRUE/FALSE}. Use this if there is more than one node.
#'
#'   \item \code{getNodeNameSingle(returnParams)}. Return the name of a
#'   single parameter/random effect node, according to \code{returnParams = 
#'   TRUE/FALSE}. Use this if there is only one node.
#'
#'   \item \code{checkInnerConvergence(message)}. Checks whether all internal
#'   optimizers converged. Returns a zero if everything converged and one
#'   otherwise. If \code{message = TRUE}, it will print more details about
#'   convergence for each conditionally independent set.
#'
#'   \item \code{gr_logLik(p, trans)}. Gradient of the (approximated)
#'   marginal log-likelihood at parameter value \code{p}. Argument \code{trans} 
#'   is similar to that in \code{calcLaplace}. If there are multiple parameters,
#'   the vector \code{p} is given in the order of parameter names returned by 
#'   \code{getNodeNamesVec(returnParams=TRUE)}.
#'
#'   \item \code{gr_Laplace(p, trans)}. This is the same as \code{gr_logLik}.
#'
#'   \item \code{otherLogLik(p)}. Calculate the \code{calcNodesOther}
#'   nodes, which returns the log-likelihood of the parts of the model that are
#'   not included in the Laplace or AGHQ approximation.
#'
#'   \item \code{gr_otherLogLik(p)}. Gradient (vector of derivatives with
#'   respect to each parameter) of \code{otherLogLik(p)}. Results should
#'   match \code{gr_otherLogLik_internal(p)} but may be more efficient after
#'   the first call.
#'
#' }
#'
#' @section Internal or development methods
#'
#' Some methods are included for calculating the (approximate) marginal log
#' posterior density by including the prior distribution of the parameters. This
#' is useful for finding the maximum a posteriori probability (MAP) estimate.
#' Currently these are provided for point calculations without estimation methods.
#'
#' \itemize{
#'
#'   \item \code{calcPrior_p(p)}. Log density of prior distribution.
#'
#'   \item \code{calcPrior_pTransformed(pTransform)}. Log density of prior distribution on transformed scale, includes the Jacobian.
#'
#'   \item \code{calcPostLogDens(p)}. Marginal log posterior density in terms of the parameter p.
#'
#'   \item \code{calcPostLogDens_pTransformed (pTransform)}. Marginal log posterior density in terms of the transformed
#'   parameter, which includes the Jacobian transformation.
#'
#'   \item \code{gr_postLogDens_pTransformed(pTransform)}. Graident of marginal log posterior density on the transformed scale. 
#'   Other available options that are used in the derivative for more flexible include \code{logDetJacobian(pTransform)} and
#'   \code{gr_logDeJacobian(pTransform)}, as well as \code{gr_prior(p)}.
#' }
#' 
#' Finally, methods that are primarily for internal use by other methods include:
#'
#' \itemize{
#'
#'    \item \code{gr_logLik_pTransformed}. Gradient of the Laplace
#'     approximation (\code{calcLogLik_pTransformed(pTransform)}) at transformed 
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
#'   \item \code{one_time_fixes()}. Users never need to run this. Is is called
#'   when necessary internally to fix dimensionality issues if there is only
#'   one parameter in the model.
#'
#'   \item \code{calcLogLik_pTransformed(pTransform)}. Laplace approximation at
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
#'   \item \code{cache_outer_logLik(logLikVal)}. Save the marginal log likelihood value 
#'   to the inner Laplace mariginlization functions to track the outer maximum internally.
#'
#'   \item \code{reset_outer_inner_logLik()}. Reset the internal saved maximum marginal log likelihood.
#'
#'   \item \code{get_inner_cholesky(atOuterMode = integer(0, default = 0))}. Returns the cholesky
#'   of the negative Hessian with respect to the random effects. If \code{atOuterMode = 1} then returns
#'   the value at the overall best marginal likelihood value, otherwise \code{atOuterMode = 0} returns the last.
#' 
#'   \item \code{get_inner_mode(atOuterMode = integer(0, default = 0))}. Returns the mode of the random effects
#'   for either the last call to the innner quadrature functions (\code{atOuterMode = 0} ), or the last best
#'   value for the marginal log likelihood, \code{atOuterMode = 1}.
#'
#' }
#'
#' @author Wei Zhang, Perry de Valpine, Paul van Dam-Bates
#' 
#' @name laplace
#' 
#' @aliases Laplace buildLaplace AGHQuad buildAGHQuad AGHQ
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
#' allres <- CpumpLaplace$summary(MLEres, randomEffectsStdError = TRUE)
#'
#' # Change the settings and also illustrate runLaplace
#' CpumpLaplace$updateSettings(innerOptimMethod = "nlminb", outerOptimMethod = "nlminb")
#' newres <- runLaplace(CpumpLaplace)
#'
#' # Illustrate use of the component log likelihood and gradient functions to
#' # run an optimizer manually from R.
#' # Use nlminb to find MLEs
#' MLEres.manual <- nlminb(c(0.1, 0.1),
#'                         function(x) -CpumpLaplace$calcLogLik(x),
#'                         function(x) -CpumpLaplace$gr_Laplace(x))
#' }
#'
#' @references
#'
#' Kass, R. and Steffey, D. (1989). Approximate Bayesian inference in
#' conditionally independent hierarchical models (parametric empirical Bayes
#' models). \emph{Journal of the American Statistical Association}, 84(407),
#' 717-726.
#' 
#' Liu, Q. and Pierce, D. A. (1994). A Note on Gauss-Hermite Quadrature. \emph{Biometrika}, 81(3) 624-629.
#'
#' Jackel, P. (2005). A note on multivariate Gauss-Hermite quadrature. London: \emph{ABN-Amro. Re.}
#'
#' Skaug, H. and Fournier, D. (2006). Automatic approximation of the marginal
#' likelihood in non-Gaussian hierarchical models. \emph{Computational
#' Statistics & Data Analysis}, 56, 699-709.
#'
NULL

# The following code takes as input a compiled Laplace approximation and returns
# a list of functions sharing an environment that provide access to the pieces
# of Laplace approximation. The main trick is to build callable interfaces to
# the AGHQuad_nfl elemensts (conditionally independent Laplace approx's), which
# are nested and so not interfaced by default. The original purpose was to
# experiment with inner (and outer) optimization methods. This is very useful
# but is not a package feature, so I am leaving the source code on display here
# but commenting it out.
## laplaceRpieces <- function(cLaplace) {
##   # limited to methodID==2
##   # uses -logLik as the working sign, so
##   # optimization will be minimization instead of maximization
##   cLaplace$one_time_fixes()
##   RLaplace <- cLaplace$Robject
##   cModel <- RLaplace$model$CobjectInterface
##   paramNodes <- RLaplace$paramNodes
##   param_values <- function(v) {
##     if(missing(v)) return(values(cModel, paramNodes))
##     else values(cModel, paramNodes) <- v
##   }
##   promoteCallable <- function(RoneLaplace, modify=TRUE) {
##     # This function is modified from the cppDef for nimbleFunctions
##     # where it is very rarely used (see comment there).
##     # There is a bug there because the indexing of existingExtPtrs is not set up
##     # so I modify here to fix that. In future, we could fix this small bug,
##     # either in promoteCallable or in the multi interface getExtPtrs method
##     # N.B. By default this modifies its argument by updating
##     # its .CobjectInterface
##     RCobj <- nimble:::nf_getRefClassObject(RoneLaplace)
##     oldCobjectInterface <- RCobj$.CobjectInterface
##     if(!is.list(oldCobjectInterface)) return(oldCobjectInterface)
##     extPtrs <- oldCobjectInterface[[1]]$getExtPtrs(oldCobjectInterface[[2]])
##     extPtrTypeIndex <- oldCobjectInterface[[1]]$extPtrTypeIndex
##     existingExtPtrs <- vector('list', length(extPtrTypeIndex))
##     existingExtPtrs[[1]] <- extPtrs[[1]]
##     existingExtPtrs[[ extPtrTypeIndex['NamedObjects']  ]] <- extPtrs[[2]]
##     thisDll <- oldCobjectInterface[[1]]$dll
##     nimbleProject <- oldCobjectInterface[[1]]$compiledNodeFun$nimbleProject
##     Rgenerator <- oldCobjectInterface[[1]]$compiledNodeFun$Rgenerator
##     newCobjectInterface <- Rgenerator(RoneLaplace, thisDll,
##                                       project = nimbleProject, existingExtPtrs = existingExtPtrs)
##     RCobj$.CobjectInterface <- newCobjectInterface
##     newCobjectInterface
##   }
##   # make sure all AGHQs (conditionally independent) are promoted
##   # to having a fully callable interface object
##   AGHQ_list_ <- RLaplace$AGHQuad_nfl$contentsList |> lapply(promoteCallable)
##   AGHQ_list_ |> lapply(\(x) x$one_time_fixes())
##   # Also make sure the parameter transformation is fully callable
##   outerParamsTransform <- promoteCallable(RLaplace$paramsTransform)
##   # reTrans refers to random effects in transformed (unconstrained) coordinates
##   # Set up some object to manage information for inner
##   # optimizations.
##   num_condIndSets <- length(AGHQ_list_)
##   # list of last optimizer results
##   last_inner_opt_list <- vector('list', num_condIndSets)
##   # list of last *best* value of inner optima
##   best_inner_opt_list <- seq_along(AGHQ_list_) |> lapply(\(x) list(value = Inf))
##   # list of constant values for reTrans to use for initializing
##   #  inner optimizations with default option and mode "constant"
##   constant_reTrans_list <- AGHQ_list_ |> lapply(
##     \(x) {
##       startID <- x$startID
##       x$startID <- 3
##       ans <- x$get_reInitTrans()
##       x$startID <- startID
##       ans
##     })
##   # list of minimum (last best) negative inner logLik value,
##   # which correspond to last_best_reTrans_list
##   best_inner_p_list <- constant_reTrans_list |>
##     lapply(\(x) rep(Inf, length(x)))
##   last_inner_p_list <- best_inner_p_list
##   reset_last_best <- function(i) {
##     if(missing(i)) i <- seq_along(AGHQ_list_)
##     for(ii in i) {
##       best_inner_opt_list[[ii]] <- list(value = Inf)
##       last_inner_opt_list[[ii]] <- list()
##       best_inner_p_list[[ii]] <- rep(Inf, length(best_inner_p_list[[ii]]))
##       last_inner_p_list[[ii]] <- best_inner_p_list[[ii]]
##     }
##   }
##   # current value of outer params
##   current_params <- numeric()
##   # current index of conditionally independent set being used
##   current_condIndSet <- 1
##   # default inner optimizer
##   default_inner_opt_fn <- \(re, fn, gr, he) {
##     nimOptim(re, fn, gr, method = "nlminb")
##   }
##   inner_opt_fn_ <- default_inner_opt_fn
##   inner_opt_fn <- function(f) {
##     if(missing(f)) return(inner_opt_fn_)
##     inner_opt_fn_ <<- f
##     f
##   }
##   # outer optimizer
##   default_outer_opt_fn <- \(p, fn, gr) {
##     optim(p, fn, gr, method = "BFGS")
##   }
##   outer_opt_fn_ <- default_outer_opt_fn
##   outer_opt_fn <- function(f) {
##     if(missing(f)) return(outer_opt_fn_)
##     outer_opt_fn_ <<- f
##     f
##   }
##   # access the list of conditionally independent AGHQs
##   AGHQ_list <- function() AGHQ_list_
##   # Objects for controlling initialization of inner optimization:
##   # Three modes are available in the default method.
##   reInitTrans_mode_ <- "constant" # or "last" or "last.best"
##   # function to set these
##   reInitTrans_mode <- function(mode) {
##     if(missing(mode)) return(reInitTrans_mode_)
##     reInitTrans_mode_ <<- mode
##   }
##   default_reInitTrans_fn <- function(AGHQobj, i) {
##     optStart <- switch(reInitTrans_mode_,
##                        last = last_inner_opt_list[[i]]$par,
##                        last.best = best_inner_opt_list[[i]],
##                        constant = constant_reTrans_list[[i]])
##     optStart
##   }
##   reInitTrans_fn_ <- default_reInitTrans_fn
##   reInitTrans_fn <- function(f) {
##     if(missing(f)) return(reInitTrans_fn_)
##     reInitTrans_fn_ <<- f
##     return(f)
##   }
##   # neg inner logLik for one conditionally independent set
##   inner_negLogLik <- function(reTrans, i = current_CondIndSet) {
##     -AGHQ_list_[[i]]$inner_logLik(reTrans)
##   }
##   # neg gradient of inner logLik for one conditionally independent set
##   gr_inner_negLogLik <- function(reTrans, i = current_CondIndSet) {
##     -AGHQ_list_[[i]]$gr_inner_logLik(reTrans)
##   }
##   # neg Hessian of inner logLik for one conditionally independent set
##   he_inner_negLogLik <- function(reTrans, i = current_CondIndSet, p = current_params) {
##     AGHQ_list_[[i]]$negHess(p, reTrans)
##   }
##   closure <- environment()
##   # minimize neg inner logLik
##   update_min_inner_negLogLik <- function(p,
##                                          reInitTrans,
##                                          i = current_CondIndSet,
##                                          inner_opt_fn,
##                                          reInitTrans_fn) {
##     optRes <- min_inner_negLogLik(p, reInitTrans, i, inner_opt_fn, reInitTrans_fn)
##     last_inner_opt_list[[i]] <- optRes
##     best_inner_p_list[[i]] <- p
##     optRes
##   }
##   min_inner_negLogLik <- function(p, reInitTrans, i = current_CondIndSet,
##                                   inner_opt_fn,
##                                   reInitTrans_fn) {
##     if(missing(inner_opt_fn)) inner_opt_fn <- get("inner_opt_fn", envir = closure)()
##     if(missing(reInitTrans_fn)) reInitTrans_fn <- get("reInitTrans_fn", envir = closure)()
##     innerObj <- AGHQ_list_[[i]]
##     if(missing(reInitTrans)) reInitTrans <- reInitTrans_fn(innerObj, i)
##     # The 1D vs nD versions different in set_params method.
##     if(length(reInitTrans) > 1) {
##       innerObj$set_params(p)
##     } else {
##       param_values(p)
##       paramDeps <- RLaplace$AGHQuad_nfl[[i]]$paramDeps
##       cModel$calculate(paramDeps)
##     }
##     fn_init <- inner_negLogLik(reInitTrans, i)
##     current_params <<- p
##     current_CondIndSet <<- i
##     if(is.nan(fn_init) || is.na(fn_init) || fn_init == Inf || fn_init == -Inf) {
##       ans <- list(par = reInitTrans, value = Inf, convergence = -1)
##       return(ans)
##     }
##     if(length(reInitTrans) > 1) {
##       if(innerObj$gr_inner_logLik_first) {
##         innerObj$gr_inner_logLik_force_update <- TRUE
##         innerObj$gr_inner_logLik(reInitTrans)
##         innerObj$gr_inner_logLik_first <- FALSE
##         innerObj$gr_inner_logLik_force_update <- FALSE
##       }
##     }
##     optRes <- inner_opt_fn(reInitTrans, fn = inner_negLogLik,
##                            gr = gr_inner_negLogLik, he = he_inner_negLogLik)
##     optRes
##   }
##   # get inner opt result
##   last_inner_opt <- function(i = current_CondIndSet) {
##     last_inner_opt_list[[i]]
##   }
##   # do one conditionally independent Laplace approx
##   one_negLaplace <- function(p, reInitTrans, i = current_CondIndset,
##                              inner_opt_fn,
##                              reInitTrans_fn,
##                              opt) {
##     if(missing(opt))
##       if(any(p!=last_inner_p_list[[i]])) {
##         if(missing(inner_opt_fn)) inner_opt_fn <- get("inner_opt_fn", envir = closure)()
##         if(missing(reInitTrans_fn)) reInitTrans_fn <- get("reInitTrans_fn", envir = closure)()
##         opt <- update_min_inner_negLogLik(p, reInitTrans, i,
##                                           inner_opt_fn, reInitTrans_fn)
##       } else {
##         opt <- last_inner_opt_list[[i]]
##       }
##     reTransform <- opt$par
##     logdetNegHessian <- AGHQ_list_[[i]]$logdetNegHess(p, reTransform)
##     nreTrans <- length(reTransform)
##     ans <- opt$value + 0.5 * logdetNegHessian - 0.5 * nreTrans * log(2*pi)
##     if(ans < best_inner_opt_list[[i]]$value) {
##       best_inner_opt_list[[i]] <- ans
##       best_inner_p_list[[i]] <- p
##     }
##     ans
##   }
##   one_gr_negLaplace <- function(p, reInitTrans, i = current_CondIndset,
##                                 inner_opt_fn,
##                                 reInitTrans_fn,
##                                 opt) {
##     if(missing(opt))
##       if(any(p!=last_inner_p_list[[i]])) {
##         if(missing(inner_opt_fn)) inner_opt_fn <- get("inner_opt_fn", envir = closure)()
##         if(missing(reInitTrans_fn)) reInitTrans_fn <- get("reInitTrans_fn", envir = closure)()
##         opt <- update_min_inner_negLogLik(p, reInitTrans, i,
##                                           inner_opt_fn, reInitTrans_fn)
##       } else {
##         opt <- last_inner_opt_list[[i]]
##       }
##     innerObj <- AGHQ_list_[[i]]
##     reTransform <- opt$par
##     negHessian <- innerObj$negHess(p, reTransform)
##     invNegHessian <- inverse(negHessian)
##     grlogdetNegHesswrtp <- innerObj$gr_logdetNegHess_wrt_p_internal(p, reTransform)
##     grlogdetNegHesswrtre <- innerObj$gr_logdetNegHess_wrt_re_internal(p, reTransform)
##     hesslogLikwrtpre <- innerObj$hess_joint_logLik_wrt_p_wrt_re_internal(p, reTransform)
##     ans <- -innerObj$gr_joint_logLik_wrt_p_internal(p, reTransform) +
##       0.5 * (grlogdetNegHesswrtp + (grlogdetNegHesswrtre %*% invNegHessian) %*% t(hesslogLikwrtpre))
##     ans[1,]
##   }
##   negLaplace <- function(p, trans = FALSE,
##                          reInitTrans,
##                          inner_opt_fn,
##                          reInitTrans_fn) {
##     if(missing(inner_opt_fn)) inner_opt_fn <- get("inner_opt_fn", envir = closure)()
##     if(missing(reInitTrans_fn)) reInitTrans_fn <- get("reInitTrans_fn", envir = closure)()
##     if(trans)
##       p <- outerParamsTransform$inverseTransform(p)
##     ans <- 0
##     if(cLaplace$num_calcNodesOther > 0) ans <- -cLaplace$otherLogLik(p)
##     missing_reInitTrans <- missing(reInitTrans)
##     for(i in seq_along(AGHQ_list_)) {
##       if(missing_reInitTrans) reIT <- reInitTrans_fn(AGHQlist[[i]], i)
##       else reIT <- reInitTrans[[i]]
##       one_ans <- one_negLaplace(p, reIT, i, inner_opt_fn, reInitTrans_fn)
##       ans <- ans + one_ans
##     }
##     if(is.nan(ans) || is.na(ans)) ans <- -Inf
##     ans
##   }
##   gr_negLaplace <- function(p, trans = FALSE,
##                          reInitTrans,
##                          inner_opt_fn,
##                          reInitTrans_fn,
##                          reset = reset_last_best) {
##     if(missing(inner_opt_fn)) inner_opt_fn <- get("inner_opt_fn", envir = closure)()
##     if(missing(reInitTrans_fn)) reInitTrans_fn <- get("reInitTrans_fn", envir = closure)()
##     if(trans) {
##       pDerivs <- cLaplace$derivs_pInverseTransform(p, 0:1)
##       p <- outerParamsTransform$inverseTransform(pDerivs$value)
##     }
##     if(cLaplace$num_calcNodesOther > 0) ans <- -cLaplace$gr_otherLogLik(p)
##     else ans <- rep(0, length(p))
##     missing_reInitTrans <- missing(reInitTrans)
##     for(i in seq_along(AGHQ_list_)) {
##       if(missing_reInitTrans) reIT <- reInitTrans_fn(AGHQlist[[i]], i)
##       else reIT <- reInitTrans[[i]]
##       one_ans <- one_gr_negLaplace(p, reIT, i, inner_opt_fn, reInitTrans_fn)
##       ans <- ans + one_ans
##     }
##     if(trans) {
##       ans <- (ans %*% pDerivs$jacobian)[1,]
##     }
##     ans
##   }
##   findMLE <- function(pStart,
##                       outer_opt_fn,
##                       inner_opt_fn,
##                       reInitTrans_fn,
##                       reset = reset_last_best) {
##     if(missing(outer_opt_fn)) outer_opt_fn <- get("outer_opt_fn", envir = closure)()
##     if(!missing(inner_opt_fn)) get("inner_opt_fn", envir = closure)(inner_opt_fun)
##     if(!missing(reInitTrans_fn)) get("reInitTrans_fn", envir = closure)(reInitTrans_fn)
##     if(is.function(reset)) reset()
##     if(missing(pStart))
##       pStart <- param_values()
##     pStartTransform <- outerParamsTransform$transform(pStart)
##     optRes <- outer_opt_fn(pStartTransform, \(p) negLaplace(p, TRUE), \(p) gr_negLaplace(p, TRUE))
##     if(optRes$convergence != 0)
##       warning("Warning: Bad outer convergence")
##     optRes$par <- outerParamsTransform$inverseTransform(optRes$par)
##     return(optRes)
##   }
##   list(promoteCallable = promoteCallable,
##        param_values = param_values,
##        inner_opt_fn = inner_opt_fn,
##        outer_opt_fn = outer_opt_fn,
##        AGHQ_list = AGHQ_list,
##        reInitTrans_mode = reInitTrans_mode,
##        reInitTrans_fn = reInitTrans_fn,
##        inner_negLogLik = inner_negLogLik,
##        gr_inner_negLogLik = gr_inner_negLogLik,
##        he_inner_negLogLik = he_inner_negLogLik,
##        min_inner_negLogLik = min_inner_negLogLik,
##        last_inner_opt = last_inner_opt,
##        one_negLaplace = one_negLaplace,
##        one_gr_negLaplace = one_gr_negLaplace,
##        negLaplace = negLaplace,
##        gr_negLaplace = gr_negLaplace,
##        outerParamsTransform = outerParamsTransform,
##        findMLE = findMLE
##        )
## }
