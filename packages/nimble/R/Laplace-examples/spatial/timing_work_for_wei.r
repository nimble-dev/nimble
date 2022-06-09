# This is a clean file to give to Wei with some efficiency improvements.

rm(list=ls())
library(nimble)

source("spatialTMB.R") # We use this to compare to TMB's performance
TMB_time # 0.71 on Perry's machine

# This is from Wei originally
source("spatial_data.R")
dd <- as.matrix(dist(Z))

# We recently made all of these use[X]Atomic options TRUE by default,
# but here they are to be explicit
nimbleOptions(useADmatInverseAtomic = TRUE)
nimbleOptions(useADcholAtomic = TRUE)
nimbleOptions(useADmatMultAtomic = TRUE)
nimbleOptions(useADsolveAtomic = TRUE)

# This option is what gives multi-level access such as CspLaplace$laplace_nfl[[1]]$negHess_inner_logLik.
# Otherwise by default access is built only one layer deep.
# (An eventual user would not turn on this option. It is for development.)
nimbleOptions(buildInterfacesForCompiledNestedNimbleFunctions = TRUE) # for getting times from Laplace

## Model code
spatialCode <- nimbleCode({
  ## Priors
  for(i in 1:2){
    b[i] ~ dnorm(0, sd = 100)
  }
  # The parameter transformation system should allow this to be more general.
  # Our TMB code for this example uses a instead of log_a and puts a constraint keeping a > 0
  # in the call to nlminb.
  # Here I changed to log_a so that a = exp(log_a) will be > 0,
  #  since we don't automatically set up box constraints.
  # a ~ dnorm(0, sd = 100)
  log_a ~ dnorm(0, sd = 100)
  a <- exp(log_a)
  log_sigma ~ dnorm(0, sd = 100)

  ## Covariance matrix
  for(i in 1:n){
    cov[i, i] <- 1
    for(j in 1:(i-1)){
      cov[i, j] <- exp(-a * dd[i, j])
      cov[j, i] <- cov[i, j]
    }
  }
  ## Zero mean
  for(i in 1:n){
    mean[i] <- 0
  }
  u[1:n] ~ dmnorm(mean[1:n], cov = cov[1:n, 1:n])

  ## Multivariate Poisson
  eta[1:n] <- X[,] %*% b[1:2] + exp(log_sigma) * u[1:n]
  for(i in 1:n){
    y[i] ~ dpois(exp(eta[i]))
  }
})

# This is an experiment.  This does not crash.  This is arguably the "right" one.
spatialGLMM <- nimbleModel(spatialCode,
                           name = "spatialGLMM",
                           constants = list(n = n, dd = dd, X = X),
                           data = list(y = y),
                           inits = list(b = c(1, 1), log_a = log(2), log_sigma = -1, u = rep(1, n)),
                           buildDerivs = TRUE
                           )

CspatialGLMM <- compileNimble(spatialGLMM)
CspatialGLMM$calculate() # check that it works

paramNodes <- c("log_a", "b", "log_sigma") # In same order as our TMB code, which confused me.
spatialGLMM$getNodeNames(topOnly = TRUE, stochOnly = TRUE) # Check what would be the default for paramNodes, but in different order

#source("Laplace_dev.R") # Load development version of Laplace
spLaplace <- buildLaplace(spatialGLMM, paramNodes = paramNodes)
CspLaplace <- compileNimble(spLaplace, project = spatialGLMM)

# Get some values for params.  Since model was just built, these are simply the inits.
p <- values(spatialGLMM, paramNodes)
tmbp <- c(exp(p[1]), p[2:4]) # tmb code has a where nimble code has log_a
# obj$fn(tmbp) #matches
# obj$gr(tmbp) # first element is divided by tmbp[1] relative to nimble

# Try one Laplace, which records the tape(s)
CspLaplace$Laplace2(p)   # record tape(s)
obj$fn(tmbp) #matches (up to flipped sign)
# Try one gradient, which records the tape(s), in case previous line wasn't run
CspLaplace$gr_Laplace2(p) # record tape(s)
obj$gr(tmbp) # first element is divided by tmbp[1] relative to nimble, due to exp(p[1])
# Try getting MLE
system.time(nimOpt <- CspLaplace$LaplaceMLE2(p)) # Currently there is a cat statement to show p[1:4] for each call.
# 2.3 s on Perry's machine
nimOpt$par
opt$par    # matches. exp(opt$par[1]) matches nimOpt$par[1]
nimOpt$value
opt$objective

# Look at times for inner optimizations
CspLaplace$laplace_nfl[[1]]$inner_logLik_times <- rep(0, 1000)
CspLaplace$laplace_nfl[[1]]$i_inner_time <- 1
system.time(nimOpt <- CspLaplace$LaplaceMLE2(p))
inner_logLik_times <- CspLaplace$laplace_nfl[[1]]$inner_logLik_times
num_inner_logLik_times <- CspLaplace$laplace_nfl[[1]]$i_inner_time - 1
inner_logLik_times[1:num_inner_logLik_times]
sum(inner_logLik_times[1:num_inner_logLik_times])
inner_logLik_times[6] / sum(inner_logLik_times[1:num_inner_logLik_times]) # A single inner opt is over half the cost of these 90 inner opts
sort(inner_logLik_times[1:num_inner_logLik_times])
hist(inner_logLik_times)


# How does nlminb for outer optimization do?
# Also I'll look at times for inner optimizations for this
CspLaplace$laplace_nfl[[1]]$inner_logLik_times <- rep(0, 1000)
CspLaplace$laplace_nfl[[1]]$i_inner_time <- 1
system.time(nimOpt_nlmb <- nlminb(p,
                                  function(p) -CspLaplace$Laplace2(p),
                                  function(p) -CspLaplace$gr_Laplace2(p))) # 1.5 s on Perry's machine
nimOpt$par # These match
num_inner_logLik_times <- CspLaplace$laplace_nfl[[1]]$i_inner_time - 1
inner_logLik_times[1:num_inner_logLik_times]
sum(inner_logLik_times[1:num_inner_logLik_times])
inner_logLik_times[6] / sum(inner_logLik_times[1:num_inner_logLik_times]) # Again the 6th call is very expensive
sort(inner_logLik_times[1:num_inner_logLik_times])


#########################
# Study efficiency

# The basic approach will be to time 100 calls to value, gradient, or hessian.
# Using 100 calls simply provides more averaging.
nimp <- p             # Use original p (some other choice should be fine too)
re <- spatialGLMM$u   # Use random effects, which happen to be in the model.extract

# This will focus on the inner log likelihood value, gr, and hess because those are the ones
# used for inner optimization, and that is the biggest cost of the whole Laplace MLE algorithm.
#
#####
# Joint log likelihood
# This is not relevant to inner optimization but is something I can compare to TMB
nim_one_fn_joint <- function(p, re) { ## assumes p has been set in model, which is how it is used
  -CspLaplace$laplace_nfl[[1]]$joint_logLik(p, re)
}
nim_one_fn_joint(p, re)
system.time(for(i in 1:100) nim_one_fn_joint(p, re)) # 0.033 s on Perry's machine

# Compare to TMB.  This is how I extracted TMB's joint logLik call.
# I can verify the numbers, but it is harder to know if there are any
# hidden efficiencies I'm missing by taking it out of its workflow.
TMB_one_fn <- function(p, re, obj) {
  par <- obj$env$par
  par[obj$env$random] <- re
  par[-obj$env$random] <- p
  obj$env$f(par, order = 0)
}
nimp <- p
tmbp <- c(exp(p[1]), p[2:4])
re <- spatialGLMM$u
TMB_one_fn(tmbp, re, obj)
system.time(for(i in 1:100) TMB_one_fn(p, re, obj)) # 0.055 s on Perry's machine

####
# Inner log likelihood
nim_one_fn_inner <- function(re) { ## assumes p has been set in model, which is how it is used
  -CspLaplace$laplace_nfl[[1]]$inner_logLik(re)
}
nim_one_fn_inner(re)
system.time(for(i in 1:100) nim_one_fn_inner(re)) # 0.017 s on Perry's machine



####
# Gradient of inner log likelihood
## This is what I think I figured out for TMB.
## It is a gradient wrt p and re, with only re elements used.
## That is what I saw in the TMB code.
## Again it is hard to unearth exactly the calling context or know how to treat
## potential R overhead.
TMB_one_inner_gr <- function(p, re, obj) {
  par <- obj$env$par
  par[obj$env$random] <- re
  par[-obj$env$random] <- p
  obj$env$f(par, order = 1)[obj$env$random] # equivalent of switch statement in ff: f0
}
gradTMBcheck <- TMB_one_inner_gr(tmbp, re, obj)

nim_one_gr_inner <- function(re) { ## assumes p has been set in model, which is how it is used
  -CspLaplace$laplace_nfl[[1]]$gr_inner_logLik(re)
}
CspLaplace$laplace_nfl[[1]]$set_params(p) # only needed once before calls with any re value
gradNimCheck <- nim_one_gr_inner(re)
max(abs(gradTMBcheck - gradNimCheck)) # They match

system.time(for(i in 1:100) TMB_one_inner_gr(tmbp, re, obj))  # 0.119 on Perry's machine
system.time(for(i in 1:100) nim_one_gr_inner(re))             # 0.024 on Perry's machine (faster than TMB)

#####
## Hessian of inner log likelihood
## Again this is the best I extracted from TMB
TMB_one_hess_inner <- function(p, re, obj) {
  par <- obj$env$par
  par[obj$env$random] <- re
  par[-obj$env$random] <- p
  obj$env$spHess(par, random = TRUE)
}
tmbhessCheck <- TMB_one_hess_inner(tmbp, re, obj)

nim_one_negHess_inner <- function(re) {
  ## This uses a triple-tape: derivs -> jac -> jac ==> Hessian
  ## We also tried double-tape straight to Hessian: derivs -> hess == > Hessian.  It was slower.
  CspLaplace$laplace_nfl[[1]]$negHess_inner_logLik(re) ## Already negative
}
# This is not already integrated into max_inner_logLik, and more intricate setup is needed:
CspLaplace$laplace_nfl[[1]]$set_params(nimp) # We need to set the outer parameters and calculated graph intermediates
CspLaplace$laplace_nfl[[1]]$record_negHess_inner_logLik(re) ## We need to record with updates on, then turn them off
nimhessCheck <- nimhess_inner <- nim_one_negHess_inner(re) # Now tis will run without wasting time on updates on each call
max(tmbhessCheck - nimhessCheck) # They match


system.time(for(i in 1:100) TMB_one_hess_inner(tmbp, re, obj))  # 0.056
system.time(for(i in 1:100) nim_one_negHess_inner(re))        # 0.127 on Perry's machine (slower than TMB)

#############################
## Try optimization

# Here is an outline of setting up different optimizations.

# I've noticed that some parameters give really bad performance,
# presumably because the arg max of the random effects is far away.
# So let's set up "easy" and "hard" parameters

p_easy <- p
p_hard <- c(-20.2624, 3.08485, -0.0406796, -5.48199) # Insanely bad parameters used in fifth BFGS step of LaplaceMLE2 above

# Our Laplace's own
# The initial re values will be:
reInit <- CspLaplace$laplace_nfl[[1]]$get_reInitTrans()
# These can also be set using
CspLaplace$laplace_nfl[[1]]$set_reInit(reInit)
# Note there is a parameter transformation involved but it is Identity for this case, so can be ignored.

# Nimble's Laplace's own inner max:
CspLaplace$laplace_nfl[[1]]$max_inner_logLik(p_easy) # The reInit was already at the max, so let's make it harder
CspLaplace$laplace_nfl[[1]]$set_reInit(rep(0, 100))
CspLaplace$laplace_nfl[[1]]$max_inner_logLik(p_easy) # That worked

system.time(CspLaplace$laplace_nfl[[1]]$max_inner_logLik(p_easy)) # 0.005 on Perry's machine
system.time(CspLaplace$laplace_nfl[[1]]$max_inner_logLik(p_hard)) # 0.485 on Perry's machine (wasn't expecting that to be faster)

########
# Pieces to try in other inner optimizers.
# Note that the functions above make everything use negative log likelihood.
# The value and gr need the negative, while the negHess, as the name implies,
# is already negative.

## nlminb
fn <- nim_one_fn_inner
gr <- nim_one_gr_inner
he <- nim_one_negHess_inner

CspLaplace$laplace_nfl[[1]]$set_params(p_easy) # We need to set outer parameters each time
nlminb(start = rep(0, 100), objective = fn, gradient = gr, hessian = he) # Answer matches

# Try nlminb WITH Hessian
CspLaplace$laplace_nfl[[1]]$set_params(p_easy)
system.time(nlminb(start = rep(0, 100), objective = fn, gradient = gr, hessian = he)) # 0.017 on Perry's machine: slower that above
CspLaplace$laplace_nfl[[1]]$set_params(p_hard)
system.time(nlminb(start = rep(0, 100), objective = fn, gradient = gr, hessian = he)) # 0.015 on Perry's machine: slower that above

# Try nlminb WITHOUT Hessian
CspLaplace$laplace_nfl[[1]]$set_params(p_easy)
system.time(nlminb(start = rep(0, 100), objective = fn, gradient = gr)) # 0.005 on Perry's machine,
CspLaplace$laplace_nfl[[1]]$set_params(p_hard)
system.time(nlminb(start = rep(0, 100), objective = fn, gradient = gr)) # 0.0055 on Perry's machine.

# Conclusion: nlminb is faster without Hessians in these examples.
# Note there is some R overhead in using nlminb, but I'm guessing not much.

#########
## Try Newton-Raphson from package maxLik
library(maxLik)
## maxLik::maxNR wants to maximize, not minimize, so flip signs again
one_nim_newton_max <- function(p, initRE) {
  CspLaplace$laplace_nfl[[1]]$set_params(p)
  maxNR(fn = CspLaplace$laplace_nfl[[1]]$inner_logLik,
        grad = CspLaplace$laplace_nfl[[1]]$gr_inner_logLik,
        hess = function(x) -CspLaplace$laplace_nfl[[1]]$negHess_inner_logLik(re),
        start = initRE, finalHessian = FALSE)
}

system.time(one_nim_newton_max(p_easy, rep(0, 100))) # 0.039 on Perry's machine
system.time(one_nim_newton_max(p_hard, rep(0, 100))) # 0.027 on Perry's machine

# CONCLUSION:
# Newton is a slower for easy p and faster for hard p
