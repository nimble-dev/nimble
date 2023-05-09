# Tests of Laplace approximation
source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

library(nimble)

# check internal consistency of optim method variants
check_laplace_alternative_methods <- function(cL, # compiled laplace algorithm
                                              cm, # compiled model
                                              m,  # original model (or list with values)
                                              opt, # possibly already-run LaplaceMLE result,
                                              methods = 1:3, # methods to check
                                              summ_orig, # summarized Laplace MLE result (original)
                                              summ_trans # summarized Laplace MLE result (transformed)
) {
  vars <- cm$getVarNames()
  reset <- function() {
    for(v in vars) cm[[v]] <- m[[v]]
  }
  if(missing(opt)) {
    reset()
    opt <- cL$LaplaceMLE()
  }
  if(missing(summ_orig)){
    summ_orig <- cL$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  }
  if(missing(summ_trans)){
    summ_trans <- cL$summary(opt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  }
  ref_method <- cL$get_method()
  for(method in methods) {
    if(method != ref_method) {
      reset()
      cL$set_method(method)
      opt_alt <- cL$LaplaceMLE()
      expect_equal(opt$par, opt_alt$par, tolerance = 0.01)
      expect_equal(opt$value, opt_alt$value, tolerance = 1e-7)
      summ_orig_alt <- cL$summary(opt_alt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
      summ_trans_alt <- cL$summary(opt_alt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
      expect_equal(summ_orig$params$estimate, summ_orig_alt$params$estimate, tol = 1e-5)
      expect_equal(summ_orig$random$estimate, summ_orig_alt$random$estimate, tol = 1e-5)
      expect_equal(summ_orig$params$stdError, summ_orig_alt$params$stdError, tol = 1e-5)
      expect_equal(summ_orig$random$stdError, summ_orig_alt$random$stdError, tol = 1e-5)
      expect_equal(summ_orig$vcov, summ_orig_alt$vcov, tol = 1e-5)
      expect_equal(summ_trans$params$estimate, summ_trans_alt$params$estimate, tol = 1e-5)
      expect_equal(summ_trans$random$estimate, summ_trans_alt$random$estimate, tol = 1e-5)
      expect_equal(summ_trans$params$stdError, summ_trans_alt$params$stdError, tol = 1e-5)
      expect_equal(summ_trans$random$stdError, summ_trans_alt$random$stdError, tol = 1e-5)
      expect_equal(summ_trans$vcov, summ_trans_alt$vcov, tol = 1e-5)
    }
  }
  invisible(NULL)
}

test_that("Laplace simplest 1D works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd = 2)
      a ~ dnorm(mu, sd = 3)
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$findMLE()
  expect_equal(opt$par, 4, tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y] = 9 + 4 = 13
  # Cov[a, y] = V[a] = 9 (not needed)
  # y ~ N(mu, 13)
  expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))
  # muhat = y = 4
  # ahat = (9*y+4*mu)/(9+4) = y = 4
  # Jacobian of ahat wrt mu is 4/13
  # Hessian of joint loglik wrt a: -(1/4 + 1/9)
  # Hessian of marginal loglik wrt mu is -1/13
  summ <- cmLaplace$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$randomEffects$estimates, 4, tol = 1e-5)
  # Covariance matrix 
  vcov <- matrix(c(1/(1/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(4/13, 1), ncol = 1) %*% (13) %*% t(matrix(c(4/13, 1), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-6)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$findMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace simplest 1D with a constrained parameter works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd = 2)
      a ~ dnorm(mu, sd = 3)
      mu ~ dexp(1.0)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  # V[a] = 9
  # V[y] = 9 + 4 = 13
  # Cov[a, y] = V[a] = 9 (not needed)
  # y ~ N(mu, 13)
  # muhat = y = 4
  # ahat = (9*y+4*mu)/(9+4) = y = 4
  # Jacobian of ahat wrt transformed param log(mu) is 4/13*mu = 4*mu/13 = 16/13
  # Hessian of joint loglik wrt a: -(1/4 + 1/9)
  # Hessian of marginal loglik wrt transformed param log(mu) is (y*mu - 2*mu*mu)/13 = -4^2/13
  # Variance of transformed param is 13/16
  expect_equal(opt$par, 4, tol = 1e-4)
  expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))
  expect_equal(opt$hessian[1,1], -4^2/13, tol = 1e-4)
  summ <- cmLaplace$summary(opt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$random$estimate, 4, tol = 1e-4)
  # Covariance matrix on transformed scale
  vcov_transform <- matrix(c(1/(1/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(16/13, 1), ncol = 1) %*% (13/16) %*% t(matrix(c(16/13, 1), ncol = 1))
  expect_equal(vcov_transform, summ$vcov, tol = 1e-4)
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  # Covariance matrix on original scale
  vcov <- diag(c(1,4)) %*% vcov_transform %*% diag(c(1,4))
  expect_equal(vcov, summ2$vcov, tol = 1e-4)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace simplest 1D (constrained) with multiple data works", {
  set.seed(1)
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, sd = 5)
      a ~ dexp(rate = exp(mu))
      for (i in 1:5){
        y[i] ~ dnorm(a, sd = 2)
      }
    }), data = list(y = rnorm(5, 1, 2)), inits = list(mu = 2, a = 1), 
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  summ <- cmLaplace$summary(opt, originalScale = FALSE, returnJointCovariance = TRUE)
  # Results are checked using those from TMB
  # TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() () 
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER(mu);
  #   PARAMETER(log_a);
  #   int n = y.size();
  #   Type a = exp(log_a); // Invserse transformation
  #   // Negative log-likelihood
  #   Type ans = -dexp(a, exp(mu), true);
  #   ans -= log_a; // logdet Jacobian of inverse transformation: exp
  #   for(int i = 0; i < n; i++){
  #     ans -= dnorm(y[i], a, Type(2), true);
  #   }
  #   return ans;
  # }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = 2, log_a = 0)
  # 
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="log_a", DLL="test")
  # tmbres <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  expect_equal(opt$par, 0.2895238)
  expect_equal(opt$value, -10.47905, tol = 1e-7)
  expect_equal(summ$random$estimate, -0.005608619)
  vcov <- matrix(c(1.414499, -1.628299, -1.628299, 2.741033), nrow = 2, byrow = TRUE)
  expect_equal(summ$vcov, vcov, 1e-6)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace simplest 1D (constrained) with deterministic intermediates and multiple data works", {
  set.seed(1)
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, sd = 5)
      a ~ dexp(rate = exp(0.5 * mu))
      for (i in 1:5){
        y[i] ~ dnorm(0.2 * a, sd = 2)
      }
    }), data = list(y = rnorm(5, 1, 2)), inits = list(mu = 2, a = 1), 
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  summ <- cmLaplace$summary(opt, originalScale = FALSE, returnJointCovariance = TRUE)
  # Results are checked using those from TMB
  # TMB cpp code:
  # #include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() () 
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER(mu);
  #   PARAMETER(log_a);
  #   int n = y.size();
  #   Type a = exp(log_a); // Invserse transformation
  #   // Negative log-likelihood
  #   Type ans = -dexp(a, exp(0.5 * mu), true);
  #   ans -= log_a; // logdet Jacobian of inverse transformation: exp
  #   for(int i = 0; i < n; i++){
  #     ans -= dnorm(y[i], 0.2 * a, Type(2), true);
  #   }
  #   ADREPORT(a);
  #   return ans;
  # }
  ## R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = 2, log_a = 0)
  # 
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="log_a", DLL="test")
  # tmbres <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  expect_equal(opt$par, -2.639534, 1e-7)
  expect_equal(opt$value, -10.47905, tol = 1e-7)
  expect_equal(summ$random$estimate, 1.603742, tol = 1e-7)
  vcov <- matrix(c(1.415167, -3.258191, -3.258191, 10.967784), nrow = 2, byrow = TRUE)
  expect_equal(summ$vcov, vcov, 1e-7)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace 1D with deterministic intermediates works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(0.2 * a, sd = 2)
      a ~ dnorm(0.5 * mu, sd = 3)
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$LaplaceMLE() # some warnings are ok here.
  expect_equal(opt$par, 40, tol = 1e-4) # 40 = 4 * (1/.2) * (1/.5)
  # V[a] = 9
  # V[y] = 0.2^2 * 9 + 4 = 4.36
  expect_equal(opt$value, dnorm(0.1*40, 0.1*40, sd = sqrt(4.36), log = TRUE))
  # y ~ N(0.2*0.5*mu, 4.36)
  # muhat = y/(0.2*0.5) = 40
  # ahat = (9*0.2*y + 4*0.5*mu)/(4+9*0.2^2) = 20
  # Jacobian of ahat wrt mu is 4*0.5/(4+9*0.2^2) = 0.4587156
  # Hessian of joint loglik wrt a: -(0.2^2/4 + 1/9)
  # Hessian of marginal loglik wrt param mu is -(0.2*0.5)^2/4.36 = -0.002293578
  summ <- cmLaplace$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$random$estimate, 20, tol = 1e-4)
  # Covariance matrix 
  vcov <- matrix(c(1/(0.2^2/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(0.4587156, 1), ncol = 1) %*% (1/0.002293578) %*% t(matrix(c(0.4587156, 1), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-4)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace 1D with a constrained parameter and deterministic intermediates works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(0.2 * a, sd = 2)
      a ~ dnorm(0.5 * mu, sd = 3)
      mu ~ dexp(1.0)
    }), data = list(y = 4), inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE() # some warnings are ok here.
  # V[a] = 9
  # V[y] = 0.2^2 * 9 + 4 = 4.36
  # y ~ N(0.2*0.5*mu, 4.36)
  # muhat = y/(0.2*0.5) = 40
  # ahat = (9*0.2*y + 4*0.5*mu)/(4+9*0.2^2) = 20
  # Jacobian of ahat wrt transformed param log(mu) is 4*0.5*mu/(4+9*0.2^2) = 18.34862
  # Hessian of joint loglik wrt a: -(0.2^2/4 + 1/9)
  # Hessian of marginal loglik wrt param mu is -(0.2*0.5)^2/4.36
  # Hessian of marginal loglik wrt transformed param log(mu) is (0.2*0.5*y*mu - 2*0.1^2*mu*mu)/4.36 = -3.669725
  expect_equal(opt$par, 40, tol = 1e-4) 
  expect_equal(opt$hessian[1,1], -3.669725, tol = 1e-3)
  expect_equal(opt$value, dnorm(0.1*40, 0.1*40, sd = sqrt(4.36), log = TRUE))
  summ <- cmLaplace$summary(opt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$random$estimate, 20, tol = 1e-4)
  # Covariance matrix on transformed scale
  vcov_transform <- matrix(c(1/(0.2^2/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(18.34862, 1), ncol = 1) %*% (1/3.669725) %*% t(matrix(c(18.34862, 1), ncol = 1))
  expect_equal(vcov_transform, summ$vcov, tol = 1e-3)
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  # Covariance matrix on original scale
  vcov <- diag(c(1,40)) %*% vcov_transform %*% diag(c(1,40))
  expect_equal(vcov, summ2$vcov, tol = 1e-3)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace 1D with deterministic intermediates and multiple data works", {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n)
        y[i] ~ dnorm(mu_y, sd = 2) # larger multiplier to amplify cov terms in result below
      mu_y <- 0.8*a
      a ~ dnorm(mu_a, sd = 3)
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }),
    data = list(y = c(4, 5, 6)),
    constants = list(n = 3),
    inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, 12.5, tol = 1e-4) # 12.5 = mean(y) * (1/.8) * (1/.5) where mean(y) = 5
  # V[a] = 9
  # V[y[i]] = 0.8^2 * 9 + 4 = 9.76
  # Cov[a, y[i]] = 0.8 * 9 = 7.2
  # Cov[y[i], y[j]] = 0.8^2 * 9 = 5.76
  Cov_ay1y2y3 <- matrix(nrow = 4, ncol = 4)
  Cov_ay1y2y3[1, 1:4] <- c(9, 7.2, 7.2, 7.2)
  Cov_ay1y2y3[2, 1:4] <- c(7.2, 9.76, 5.76, 5.76)
  Cov_ay1y2y3[3, 1:4] <- c(7.2, 5.76, 9.76, 5.76)
  Cov_ay1y2y3[4, 1:4] <- c(7.2, 5.76, 5.76, 9.76)
  Cov_y1y2y3 <- Cov_ay1y2y3[2:4, 2:4]
  chol_cov <- chol(Cov_y1y2y3)
  res <- dmnorm_chol(c(4, 5, 6), 0.8*0.5*12.5, cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # y[i] ~ N(0.4*mu, 9.76) 
  # mean(y) = 5
  # muhat = mean(y)/(0.8*0.5) = 12.5
  # ahat = (9*0.8*sum(y) + 4*0.5*mu)/(4+9*0.8^2*3) = 6.25
  # Jacobian of ahat wrt mu is 4*0.5/(4+9*0.8^2*3) = 0.09398496
  # Hessian of joint loglik wrt a: -(3*0.8^2/4 + 1/9)
  # Hessian of marginal loglik wrt mu: -0.02255639 (numerical, have not got AD work)
  summ <- cmLaplace$summary(opt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$random$estimate, 6.25, tol = 1e-6)
  # Covariance matrix 
  vcov <- matrix(c(1/(0.8^2*3/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(0.09398496, 1), ncol = 1) %*% (1/0.02255639) %*% t(matrix(c(0.09398496, 1), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-7)
  
  # check that
  mLaplaceCheck <- buildLaplace(model = m, paramNodes = 'mu', randomEffectsNodes = 'a')
  nim1D <-  mLaplace$laplace_nfl[[1]]
  expect_identical(nim1D$paramNodes, "mu")
  expect_identical(nim1D$paramDeps, "mu_a")
  expect_identical(nim1D$randomEffectsNodes, "a")
  expect_identical(nim1D$innerCalcNodes, c("a", "mu_y", "y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$calcNodes, c("mu_a", nim1D$innerCalcNodes))
  expect_identical(nim1D$inner_updateNodes, "mu_a")
  expect_identical(nim1D$inner_constantNodes, c("y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$joint_updateNodes, character())
  expect_identical(nim1D$joint_constantNodes, c("y[1]", "y[2]", "y[3]"))

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace 1D with a constrained parameter and deterministic intermediates and multiple data works", {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n)
        y[i] ~ dnorm(mu_y, sd = 2) 
      mu_y <- 0.8*a
      a ~ dnorm(mu_a, sd = 3)
      mu_a <- 0.5 * mu
      mu ~ dexp(1.0)
    }),
    data = list(y = c(4, 5, 6)),
    constants = list(n = 3),
    inits = list(a = -1, mu = 0),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, 12.5, tol = 1e-4)
  # V[a] = 9
  # V[y[i]] = 0.8^2 * 9 + 4 = 9.76
  # Cov[a, y[i]] = 0.8 * 9 = 7.2
  # Cov[y[i], y[j]] = 0.8^2 * 9 = 5.76
  # y[i] ~ N(0.4*mu, 9.76) 
  # mean(y) = 5
  # muhat = mean(y)/(0.8*0.5) = 12.5
  # ahat = (9*0.8*sum(y) + 4*0.5*mu)/(4+9*0.8^2*3) = 6.25
  # Jacobian of ahat wrt transformed param log(mu) is 4*0.5*mu/(4+9*0.8^2*3) = 1.174812
  # Hessian of joint loglik wrt a: -(3*0.8^2/4 + 1/9)
  # Hessian of marginal loglik wrt transformed param: -3.524436 (numerical, have not got AD work)
  Cov_ay1y2y3 <- matrix(nrow = 4, ncol = 4)
  Cov_ay1y2y3[1, 1:4] <- c(9, 7.2, 7.2, 7.2)
  Cov_ay1y2y3[2, 1:4] <- c(7.2, 9.76, 5.76, 5.76)
  Cov_ay1y2y3[3, 1:4] <- c(7.2, 5.76, 9.76, 5.76)
  Cov_ay1y2y3[4, 1:4] <- c(7.2, 5.76, 5.76, 9.76)
  Cov_y1y2y3 <- Cov_ay1y2y3[2:4, 2:4]
  chol_cov <- chol(Cov_y1y2y3)
  res <- dmnorm_chol(c(4, 5, 6), 0.8*0.5*12.5, cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # Check covariance matrix 
  summ <- cmLaplace$summary(opt, originalScale = FALSE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$random$estimate, 6.25, tol = 1e-6)
  # Covariance matrix on transformed scale
  vcov_transform <- matrix(c(1/(0.8^2*3/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(1.174812, 1), ncol = 1) %*% (1/3.524436) %*% t(matrix(c(1.174812, 1), ncol = 1))
  expect_equal(vcov_transform, summ$vcov, tol = 1e-6)
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  # Covariance matrix on original scale
  vcov <- diag(c(1,12.5)) %*% vcov_transform %*% diag(c(1,12.5))
  expect_equal(vcov, summ2$vcov, tol = 1e-5)
  # check that
  mLaplaceCheck <- buildLaplace(model = m, paramNodes = 'mu', randomEffectsNodes = 'a')
  nim1D <-  mLaplace$laplace_nfl[[1]]
  expect_identical(nim1D$paramNodes, "mu")
  expect_identical(nim1D$paramDeps, "mu_a")
  expect_identical(nim1D$randomEffectsNodes, "a")
  expect_identical(nim1D$innerCalcNodes, c("a", "mu_y", "y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$calcNodes, c("mu_a", nim1D$innerCalcNodes))
  expect_identical(nim1D$inner_updateNodes, "mu_a")
  expect_identical(nim1D$inner_constantNodes, c("y[1]", "y[2]", "y[3]"))
  expect_identical(nim1D$joint_updateNodes, character())
  expect_identical(nim1D$joint_constantNodes, c("y[1]", "y[2]", "y[3]"))
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace simplest 2x1D works, with multiple data for each", {
  set.seed(1)
  y <- matrix(rnorm(6, 4, 5), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) {
        mu_y[i] <- 0.8*a[i]
        for(j in 1:3)
          y[i, j] ~ dnorm(mu_y[i], sd = 2)
        a[i] ~ dnorm(mu_a, sd = 3)
      }
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }), data = list(y = y), inits = list(a = c(-2, -1), mu = 0),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, mean(y)/(0.8*0.5), tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y[i]] = 0.8^2 * 9 + 4 = 9.76
  # Cov[a, y[i]] = 0.8 * 9 = 7.2
  # Cov[y[i], y[j]] = 0.8^2 * 9 = 5.76, within a group
  Cov_A_Y <- matrix(nrow = 8, ncol = 8)
  Cov_A_Y[1, 1:8] <- c(  9,    0,  7.2,  7.2,  7.2,    0,    0,    0)
  Cov_A_Y[2, 1:8] <- c(  0,    9,    0,    0,    0,  7.2,  7.2,  7.2)
  Cov_A_Y[3, 1:8] <- c(7.2,    0, 9.76, 5.76, 5.76,    0,    0,    0)
  Cov_A_Y[4, 1:8] <- c(7.2,    0, 5.76, 9.76, 5.76,    0,    0,    0)
  Cov_A_Y[5, 1:8] <- c(7.2,    0, 5.76, 5.76, 9.76,    0,    0,    0)
  Cov_A_Y[6, 1:8] <- c(  0,  7.2,    0,    0,    0, 9.76, 5.76, 5.76)
  Cov_A_Y[7, 1:8] <- c(  0,  7.2,    0,    0,    0, 5.76, 9.76, 5.76)
  Cov_A_Y[8, 1:8] <- c(  0,  7.2,    0,    0,    0, 5.76, 5.76, 9.76)
  Cov_Y <- Cov_A_Y[3:8, 3:8]
  chol_cov <- chol(Cov_Y)
  res <- dmnorm_chol(as.numeric(t(y)), mean(y), cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # muhat = mean(y)/(0.8*0.5)
  # ahat[1] = (9*0.8*sum(y[1,]) + 4*0.5*mu)/(4+9*0.8^2*3)
  # ahat[2] = (9*0.8*sum(y[2,]) + 4*0.5*mu)/(4+9*0.8^2*3)
  # Jacobian of ahat[i] wrt mu is 4*0.5/(4+9*0.8^2*3) = 0.09398496
  # Hessian of joint loglik wrt a[i]a[i]: -(3*0.8^2/4 + 1/9); wrt a[i]a[j]: 0
  # Hessian of marginal loglik wrt mu: -0.04511278 (numerical, have not got AD work)
  muhat <- mean(y)/(0.8*0.5)
  ahat <- c((9*0.8*sum(y[1,]) + 4*0.5*muhat)/(4+9*0.8^2*3), (9*0.8*sum(y[2,]) + 4*0.5*muhat)/(4+9*0.8^2*3))
  summ <- cmLaplace$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$random$estimate, ahat, tol = 1e-6)
  # Covariance matrix 
  vcov <- diag(c(rep(1/(3*0.8^2/4 + 1/9), 2), 0)) + matrix(c(rep(0.09398496, 2), 1), ncol = 1) %*% (1/0.04511278) %*% t(matrix(c(rep(0.09398496, 2), 1), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-7)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace with 2x1D random effects needing joint integration works, without intermediate nodes", {
  set.seed(1)
  y <- matrix(rnorm(6, 4, 5), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) {
        a[i] ~ dnorm(mu_a, sd = 3)
      }
      for(j in 1:3) # Note this is different than above.
        # These are 3 observations each of 2D
        y[1:2, j] ~ dmnorm(a[1:2], cov = cov_y[1:2, 1:2])
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }),
    data = list(y = y),
    inits = list(a = c(-2, -1), mu = 0),
    constants = list(cov_y = matrix(c(2, 1.5, 1.5, 2), nrow = 2)),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, mean(y)/(0.5), tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y[1:2, i]] =  diag(2)*9 + cov_y = [ (9 + 2), 0 + 1.5; 0+1.5, (9+2)]
  # Cov[a[1], y[1,j]] = 9
  # Cov[a[1], y[2,j]] = 0
  # Cov[a[2], y[1,j]] = 0
  # Cov]a[2], y[2,j]] = 9
  # Cov[y[1,i], y[1,j]] = 9
  # Cov[y[2,i], y[2,j]] = 9
  Cov_A_Y <- matrix(nrow = 8, ncol = 8)
  Cov_A_Y[1, 1:8] <- c(  9,    0,    9,    0,    9,    0,    9,    0)
  Cov_A_Y[2, 1:8] <- c(  0,    9,    0,    9,    0,    9,    0,    9)
  Cov_A_Y[3, 1:8] <- c(  9,    0,   11,  1.5,    9,    0,    9,    0)
  Cov_A_Y[4, 1:8] <- c(  0,    9,  1.5,   11,    0,    9,    0,    9)
  Cov_A_Y[5, 1:8] <- c(  9,    0,    9,    0,   11,  1.5,    9,    0)
  Cov_A_Y[6, 1:8] <- c(  0,    9,    0,    9,  1.5,   11,    0,    9)
  Cov_A_Y[7, 1:8] <- c(  9,    0,    9,    0,    9,    0,   11,  1.5)
  Cov_A_Y[8, 1:8] <- c(  0,    9,    0,    9,    0,    9,  1.5,   11)
  Cov_Y <- Cov_A_Y[3:8, 3:8]
  chol_cov <- chol(Cov_Y)
  res <- dmnorm_chol(as.numeric(y), mean(y), cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, returnJointCovariance = TRUE)
  # Covariance matrix from TMB:
  # TMB cpp code (test.cpp) below:
  # include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() () 
  # {
  #   DATA_MATRIX(y);
  #   DATA_MATRIX(Sigma);
  #   PARAMETER(mu);
  #   PARAMETER_VECTOR(a);
  #   int i;
  #   Type ans = 0.0;
  #   // Negative log-likelihood
  #   for(i = 0; i < 2; i++){
  #     ans -= dnorm(a[i], 0.5*mu, Type(3.0), true);
  #   }
  #   vector<Type> residual(2);
  #   using namespace density;
  #   MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  #   for(i = 0; i < 3; i++)
  #   {
  #     residual = vector<Type>(y.col(i)) - a;
  #     ans += neg_log_dmvnorm(residual);
  #   }
  #   return ans;
  #   }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp") 
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y, Sigma = m$cov_y)
  # parameters <- list(mu = 0, a = c(-2, -1))
  # 
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  tmbvcov <- matrix(nrow = 3, ncol = 3)
  tmbvcov[1,] <- c(20.333333, 1.1666667, 1.1666667)
  tmbvcov[2,] <- c(1.166667, 0.6651515, 0.5015152)
  tmbvcov[3,] <- c(1.166667, 0.5015152, 0.6651515)
  ## TMB and NIMBLE have different orders of random effects and parameters
  expect_equal(summ$vcov[c(3,1,2), c(3,1,2)], tmbvcov, tol = 1e-6)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace with 2x1D random effects needing joint integration works, with intermediate nodes", {
  set.seed(1)
  y <- matrix(rnorm(6, 4, 5), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) {
        mu_y[i] <- 0.8*a[i]
        a[i] ~ dnorm(mu_a, sd = 3)
      }
      for(j in 1:3)
        y[1:2, j] ~ dmnorm(mu_y[1:2], cov = cov_y[1:2, 1:2])
      mu_a <- 0.5 * mu
      mu ~ dnorm(0, sd = 5)
    }),
    data = list(y = y),
    inits = list(a = c(-2, -1), mu = 0),
    constants = list(cov_y = matrix(c(2, 1.5, 1.5, 2), nrow = 2)),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, mean(y)/(0.8*0.5), tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y[1:2, i]] =  diag(2)*0.8^2 * 9 + cov_y = [ (5.76 + 2), 0 + 1.5; 0+1.5, (5.76+2)]
  # Cov[a[1], y[1,j]] = 0.8*9 = 7.2
  # Cov[a[1], y[2,j]] = 0
  # Cov[a[2], y[1,j]] = 0
  # Cov]a[2], y[2,j]] = 0.8*9
  # Cov[y[1,i], y[1,j]] = 0.8^2*9 = 5.76
  # Cov[y[2,i], y[2,j]] = 9
  Cov_A_Y <- matrix(nrow = 8, ncol = 8)
  Cov_A_Y[1, 1:8] <- c(  9,    0,  7.2,    0,  7.2,    0,  7.2,    0)
  Cov_A_Y[2, 1:8] <- c(  0,    9,    0,  7.2,    0,  7.2,    0,  7.2)
  Cov_A_Y[3, 1:8] <- c(7.2,    0, 7.76,  1.5, 5.76,    0, 5.76,    0)
  Cov_A_Y[4, 1:8] <- c(  0,  7.2,  1.5, 7.76,    0, 5.76,    0, 5.76)
  Cov_A_Y[5, 1:8] <- c(7.2,    0, 5.76,    0, 7.76,  1.5, 5.76,    0)
  Cov_A_Y[6, 1:8] <- c(  0,  7.2,    0, 5.76,  1.5, 7.76,    0, 5.76)
  Cov_A_Y[7, 1:8] <- c(7.2,    0, 5.76,    0, 5.76,    0, 7.76,  1.5)
  Cov_A_Y[8, 1:8] <- c(  0,  7.2,    0, 5.76,    0, 5.76,  1.5, 7.76)
  Cov_Y <- Cov_A_Y[3:8, 3:8]
  chol_cov <- chol(Cov_Y)
  res <- dmnorm_chol(as.numeric(y), mean(y), cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  expect_equal(opt$value, res)
  
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, returnJointCovariance = TRUE)
  # Covariance matrix from TMB:
  # TMB cpp code (test.cpp) below:
  # include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() () 
  # {
  #   DATA_MATRIX(y);
  #   DATA_MATRIX(Sigma);
  #   PARAMETER(mu);
  #   PARAMETER_VECTOR(a);
  #   int i;
  #   Type ans = 0.0;
  #   // Negative log-likelihood
  #   for(i = 0; i < 2; i++){
  #     ans -= dnorm(a[i], 0.5*mu, Type(3.0), true);
  #   }
  #   vector<Type> residual(2);
  #   using namespace density;
  #   MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  #   for(i = 0; i < 3; i++)
  #   {
  #     residual = vector<Type>(y.col(i)) - 0.8 * a;
  #     ans += neg_log_dmvnorm(residual);
  #   }
  #   return ans;
  #   }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp") 
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y, Sigma = m$cov_y)
  # parameters <- list(mu = 0, a = c(-2, -1))
  # 
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  tmbvcov <- matrix(nrow = 3, ncol = 3)
  tmbvcov[1,] <- c(21.645833, 1.8229167, 1.8229167)
  tmbvcov[2,] <- c(1.822917, 1.0380050, 0.7849117)
  tmbvcov[3,] <- c(1.822917, 0.7849117, 1.0380050)
  ## TMB and NIMBLE have different orders of random effects and parameters
  expect_equal(summ$vcov[c(3,1,2), c(3,1,2)], tmbvcov, tol = 1e-6)

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace with 2x2D random effects for 1D data that are separable works, with intermediate nodes", {
  set.seed(1)
  # y[i, j] is jth datum from ith group
  y <- array(rnorm(8, 6, 5), dim = c(2, 2, 2)) 
  cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
      mu_a[1] <- 0.8 * mu[1]
      mu_a[2] <- 0.2 * mu[2]
      for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
      for(i in 1:2) {
        for(j in 1:2) {
          y[1, j, i] ~ dnorm( 0.5 * a[i, 1], sd = 1.8) # this ordering makes it easier below
          y[2, j, i] ~ dnorm( 0.1 * a[i, 2], sd = 1.2)
        }
      }
    }),
    data = list(y = y),
    inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2), mu = c(0, .5)),
    constants = list(cov_a = cov_a),
    buildDerivs = TRUE
  )

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$LaplaceMLE()
  
  ## Wei: I tested this using TMB instead of the code below
  # TMB cpp code:
  # #include <TMB.hpp>
  # template<class Type>
  # Type objective_function<Type>::operator() () 
  # {
  #   DATA_ARRAY(y);
  #   DATA_MATRIX(Sigma);
  #   PARAMETER_VECTOR(mu);
  #   PARAMETER_MATRIX(a);
  #   int i, j;
  #   Type ans = 0.0;
  #   vector<Type> mu_a(2);
  #   mu_a(0) = 0.8 * mu(0);
  #   mu_a(1) = 0.2 * mu(1);
  #   // Negative log-likelihood
  #   vector<Type> residual(2);
  #   using namespace density;
  #   MVNORM_t<Type> neg_log_dmvnorm(Sigma);
  #   for(i = 0; i < 2; i++)
  #   {
  #     residual = vector<Type>(a.row(i)) - mu_a;
  #     ans += neg_log_dmvnorm(residual);
  #   }
  #   for(i = 0; i < 2; i++){
  #     for(j = 0; j < 2; j++){
  #       ans -= dnorm(y(0, j, i), 0.5*a(i, 0), Type(1.8), true);
  #       ans -= dnorm(y(1, j, i), 0.1*a(i, 1), Type(1.2), true);
  #     }
  #   }
  #   return ans;
  # }
  # TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y, Sigma = m$cov_a)
  # parameters <- list(mu = m$mu, a = m$a)
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  expect_equal(opt$par, c(12.98392, 406.04878), tol = 1e-4)
  expect_equal(opt$value, -41.86976, tol = 1e-6)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, returnJointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 6, ncol = 6)
  tmbvcov[1,] <- c(6.625000e+00, 4.687500e+00,  4.050000e+00,  4.050000e+00, -2.693817e-11, -2.695275e-11)
  tmbvcov[2,] <- c(4.687500e+00, 9.250000e+02,  2.965628e-11,  2.967848e-11,  1.800000e+02,  1.800000e+02)
  tmbvcov[3,] <- c(4.050000e+00, 2.951367e-11,  3.995242e+00,  2.484758e+00,  5.596302e-01, -5.596302e-01)
  tmbvcov[4,] <- c(4.050000e+00, 2.951367e-11,  2.484758e+00,  3.995242e+00, -5.596302e-01,  5.596302e-01)
  tmbvcov[5,] <- c(-2.691772e-11, 1.800000e+02,  5.596302e-01, -5.596302e-01,  3.684693e+01,  3.515307e+01)
  tmbvcov[6,] <- c(-2.691772e-11, 1.800000e+02, -5.596302e-01,  5.596302e-01,  3.515307e+01,  3.684693e+01)
  
  expect_equal(summ$vcov[c(5,6,1,3,2,4), c(5,6,1,3,2,4)], tmbvcov, tol = 1e-4)
  
  # For this case, we build up the correct answer more formulaically
  # Define A as the vector a[1, 1], a[1, 2], a[2, 1], a[2, 2]
  # cov_A <- matrix(0, nrow = 4, ncol = 4)
  # cov_A[1:2, 1:2] <- cov_a
  # cov_A[3:4, 3:4] <- cov_a
  # # Define Y as the vector y[1,1,1],y[2,1,1],y[1,2,1],y[2,2,1], then same with last index 2
  # # Define E[Y] as IA %*% A, where:
  # IA <- matrix(0, nrow = 8, ncol = 4)
  # IA[c(1, 3), 1] <- 0.5
  # IA[c(2, 4), 2] <- 0.1
  # IA[c(5, 7), 3] <- 0.5
  # IA[c(6, 8), 4] <- 0.1
  # 
  # # define cov_y_given_a as the Cov[Y | A]
  # cov_y_given_a <- matrix(0, nrow = 8, ncol = 8)
  # diag(cov_y_given_a) <- rep(c(1.8^2, 1.2^2), 4)
  # # And finally get cov_Y, the marginal (over A) covariance of Y
  # cov_Y <- IA %*% cov_A %*% t(IA) + cov_y_given_a
  # chol_cov <- chol(cov_Y)
  # 
  # # make a log likelihood function
  # nlogL <- function(mu) {
  #   mean_Y <- rep(c(0.8*0.5*mu[1], 0.2*0.1*mu[2]), 4)
  #   -dmnorm_chol(as.numeric(y), mean_Y, cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  # }
  # # maximize it
  # opt_manual <- optim(c(20, 100), nlogL, method = "BFGS")
  # expect_equal(opt$par, opt_manual$par, tol = 1e-4)
  # expect_equal(opt$value, -opt_manual$value, tol = 1e-5)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace with 2x2D random effects for 2D data that need joint integration works, with intermediate nodes", {
  set.seed(1)
  cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  cov_y <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  y <- rmnorm_chol(1, c(1, 1), chol(cov_y), prec_param = FALSE)
  y <- rbind(y, rmnorm_chol(1, c(1, 1), chol(cov_y), prec_param = FALSE))
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
      mu_a[1] <- 0.8 * mu[1]
      mu_a[2] <- 0.2 * mu[2]
      for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
      mu_y[1:2] <- 0.5*a[1, 1:2] + 0.1*a[2, 1:2]
      for(i in 1:2) {
        y[i, 1:2] ~ dmnorm(mu_y[1:2], cov = cov_y[1:2, 1:2])
      }
    }),
    data = list(y = y),
    inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2), mu = c(0, 0.5)),
    constants = list(cov_a = cov_a, cov_y = cov_y),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  ## Check using TMB results
  expect_equal(opt$par, c(0.5603309, 11.7064674 ), tol = 1e-4)
  expect_equal(opt$value, -4.503796, tol = 1e-7)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, returnJointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 6, ncol = 6)
  tmbvcov[1,] <- c(4.4270833,  11.111111, 1.4583333, 3.1250000, 0.6597222,  1.9097222)
  tmbvcov[2,] <- c(11.1111111, 70.833333, 2.6388889, 7.6388889, 5.8333333, 12.5000000)
  tmbvcov[3,] <- c(1.4583333,   2.638889, 1.5000000, 0.8333333, 0.7777778,  0.2777778)
  tmbvcov[4,] <- c(3.1250000,   7.638889, 0.8333333, 4.1666667, 0.2777778,  2.7777778)
  tmbvcov[5,] <- c(0.6597222,   5.833333, 0.7777778, 0.2777778, 1.5000000,  0.8333333)
  tmbvcov[6,] <- c(1.9097222,  12.500000, 0.2777778, 2.7777778, 0.8333333,  4.1666667)
  
  expect_equal(summ$vcov[c(5,6,1,3,2,4), c(5,6,1,3,2,4)], tmbvcov, tol = 1e-4)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
  ## TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() () 
  # {
  #   DATA_MATRIX(y);
  #   DATA_MATRIX(cov_a);
  #   DATA_MATRIX(cov_y);
  #   PARAMETER_VECTOR(mu);
  #   PARAMETER_MATRIX(a);
  #   int i;
  #   Type ans = 0.0;
  #   
  #   using namespace density;
  #   // Negative log-likelihood of mv normal
  #   vector<Type> mu_a(2);
  #   mu_a(0) = 0.8 * mu(0);
  #   mu_a(1) = 0.2 * mu(1);
  #   vector<Type> residual_a(2);
  #   MVNORM_t<Type> dmvnorm_a(cov_a);
  #   for(i = 0; i < 2; i++)
  #   {
  #     residual_a = vector<Type>(a.row(i)) - mu_a;
  #     ans += dmvnorm_a(residual_a);
  #   }
  #   vector<Type> mu_y(2);
  #   mu_y(0) = 0.5*a(0, 0) + 0.1*a(1, 0);
  #   mu_y(1) = 0.5*a(0, 1) + 0.1*a(1, 1);
  #   vector<Type> residual_y(2);
  #   MVNORM_t<Type> dmvnorm_y(cov_y);
  #   for(i = 0; i < 2; i++){
  #     residual_y = vector<Type>(y.row(i)) - mu_y;
  #     ans += dmvnorm_y(residual_y);
  #   }
  #   return ans;
  # }
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y,  cov_a = m$cov_a, cov_y = m$cov_y)
  # parameters <- list(mu = m$mu, a = m$a)
  # 
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
})

test_that("simple LME case works", {
  set.seed(1)
  g <- rep(1:10, each = 5)
  n <- length(g)
  x <- runif(n)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n) {
        y[i] ~ dnorm((fixed_int + random_int[g[i]]) + (fixed_slope + random_slope[g[i]])*x[i], sd = sigma_res)
      }
      for(i in 1:ng) {
        random_int[i] ~ dnorm(0, sd = sigma_int)
        random_slope[i] ~ dnorm(0, sd = sigma_slope)
      }
      sigma_int ~ dunif(0, 10)
      sigma_slope ~ dunif(0, 10)
      sigma_res ~ dunif(0, 10)
      fixed_int ~ dnorm(0, sd = 100)
      fixed_slope ~ dnorm(0, sd = 100)
    }),
    constants = list(g = g, ng = max(g), n = n, x = x),
    buildDerivs = TRUE
  )
  params <- c("fixed_int", "fixed_slope", "sigma_int", "sigma_slope", "sigma_res")
  values(m, params) <- c(10, 0.5, 3, .25, 0.2)
  m$simulate(m$getDependencies(params, self = FALSE))
  m$setData('y')
  y <- m$y
  library(lme4)
  manual_fit <- lmer(y ~ x + (1 + x || g), REML = FALSE)

  mLaplace <- buildLaplace(model = m)
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  opt <- cmLaplace$LaplaceMLE()
  nimres <- cmLaplace$summary(opt, calcRandomEffectsStdError = TRUE)
  lme4res <- summary(manual_fit)
  expect_equal(nimres$params$estimate[4:5], as.vector(lme4res$coefficients[,"Estimate"]), tol=1e-6)
  expect_equal(nimres$params$estimate[1:3], as.data.frame(VarCorr(manual_fit))[,"sdcor"], tol = 1e-5)
  expect_equal(nimres$params$stdError[4:5], as.vector(lme4res$coefficients[,"Std. Error"]), tol=1e-4)
  expect_equal(nimres$random$estimate, as.vector(t(ranef(manual_fit)$g)), tol = 1e-5)
})

test_that("simple LME with correlated intercept and slope works", {
  set.seed(1)
  g <- rep(1:10, each = 10)
  n <- length(g)
  x <- runif(n)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:n) {
        y[i] ~ dnorm((fixed_int + random_int_slope[g[i], 1]) + (fixed_slope + random_int_slope[g[i], 2])*x[i], sd = sigma_res)
      }
      cov[1, 1] <- sigma_int^2
      cov[2, 2] <- sigma_slope^2
      cov[1, 2] <- rho * sigma_int * sigma_slope
      cov[2, 1] <- rho * sigma_int * sigma_slope
      for(i in 1:ng) {
        random_int_slope[i, 1:2] ~ dmnorm(zeros[1:2], cov = cov[1:2, 1:2])
      }
      sigma_int ~ dunif(0, 10)
      sigma_slope ~ dunif(0, 10)
      sigma_res ~ dunif(0, 10)
      fixed_int ~ dnorm(0, sd = 100)
      fixed_slope ~ dnorm(0, sd = 100)
      rho ~ dunif(-1, 1)
    }),
    constants = list(g = g, ng = max(g), n = n, x = x, zeros = rep(0, 2)),
    buildDerivs = TRUE
  )
  params <- c("fixed_int", "fixed_slope", "sigma_int", "sigma_slope", "sigma_res", "rho")
  values(m, params) <- c(10, 0.5, 3, 0.25, 0.2, 0.45)
  m$simulate(m$getDependencies(params, self = FALSE))
  m$setData('y')
  y <- m$y
  library(lme4)
  manual_fit <- lmer(y ~ x + (1 + x | g), REML = FALSE)
  mLaplace <- buildLaplace(model = m)
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  opt <- cmLaplace$LaplaceMLE()
  nimres <- cmLaplace$summary(opt, calcRandomEffectsStdError = TRUE)
  lme4res <- summary(manual_fit)
  expect_equal(nimres$params$estimate[4:5], as.vector(lme4res$coefficients[,"Estimate"]), tol=1e-4)
  sdparams <- nimres$params$estimate[-c(4,5)]
  expect_equal(sdparams[c(1,2,4,3)], as.data.frame(VarCorr(manual_fit))[,"sdcor"], tol = 1e-3)
  expect_equal(nimres$params$stdError[4:5], as.vector(lme4res$coefficients[,"Std. Error"]), tol=5e-3)
  expect_equal(nimres$random$estimate, as.vector(t(ranef(manual_fit)$g)), tol = 5e-3)
})

test_that("Laplace with non-empty calcNodesNoLaplace works", {
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:3) {
        mu[i] ~ dnorm(0, sd = 10)
      }
      mu_a[1] <- mu[1] + mu[2]
      mu_a[2] <- mu[2] + mu[3]
      a[1] ~ dnorm(mu_a[1], sd = 2)
      y[1] ~ dnorm(a[1], sd = 3)
      a[2] ~ dnorm(mu_a[2], sd = 2)
      y[2] ~ dnorm(a[2], sd =3)
      y[3] ~ dnorm(mu[3], sd = 3)
    }),
    data = list(y = c(2, 3, 5)),
    inits = list(a = c(1, 2), mu = c(1, 2, 3)),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, c(4, -2, 5), tol = 1e-3)
  expect_equal(opt$value, -6.420377, tol = 1e-6)
  
  ## Check covariance matrix
  summ <- cmLaplace$summary(opt, returnJointCovariance = TRUE)
  ## TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() () 
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER_VECTOR(mu);
  #   PARAMETER_VECTOR(a);
  #   int i;
  #   // Negative log-likelihood
  #   Type ans = -dnorm(a[0], mu[0]+mu[1], Type(2.0), true);
  #   ans -= dnorm(a[1], mu[1]+mu[2], Type(2.0), true);
  #   for(i = 0; i < 2; i++){
  #     ans -= dnorm(y[i], a[i], Type(3.0), true);
  #   }
  #   ans -= dnorm(y[2], mu[2], Type(3.0), true);
  #   return ans;
  # }
  ## TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = c(1, 2, 3), a = c(1, 2))
  # 
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="a", DLL="test")
  # tmbres <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
  
  ## Covariance matrix from TMB
  tmbvcov <- matrix(nrow = 5, ncol = 5)
  tmbvcov[1,] <- c( 35, -2.20000e+01,  9.000000e+00,  9.000000e+00, -9.000000e+00)
  tmbvcov[2,] <- c(-22,  2.20000e+01, -9.000000e+00,  8.463230e-13,  9.000000e+00)
  tmbvcov[3,] <- c( 9,  -9.00000e+00,  9.000000e+00, -3.462231e-13,  3.462231e-13)
  tmbvcov[4,] <- c( 9,   8.46323e-13, -3.462231e-13,  9.000000e+00,  3.462231e-13)
  tmbvcov[5,] <- c(-9,   9.00000e+00,  3.462231e-13,  3.462231e-13,  9.000000e+00)
  
  expect_equal(summ$vcov[c(3:5, 1:2), c(3:5, 1:2)], tmbvcov, tol=1e-5)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() 
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

test_that("Laplace with 2x1D parameters (one needs transformation) and non-normal data works", {
  m <- nimbleModel(
    nimbleCode({
      mu ~ dnorm(0, sd = 10.0)
      sigma ~ dunif(0, 100)
      for (i in 1:5){
        theta[i] ~ dnorm(mu, sd = sigma)
        logit(p[i]) <- theta[i]
        y[i] ~ dbinom(10, prob = p[i])
      }
    }),
    data = list(y = c(8, 6, 5, 3, 7)),
    inits = list(mu = 1, sigma = 1, theta = rep(0, 5)),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  ## Compare with results from TMB
  expect_equal(opt$par, c(0.330241, 0.3059177), tol = 1e-4)
  expect_equal(opt$value, -9.703857, tol = 1e-6)
  ## Check covariance matrix on the transformed scale
  summ <- cmLaplace$summary(opt, originalScale = FALSE, returnJointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 7, ncol = 7)
  tmbvcov[1,] <- c(0.10337427,  0.04574391,  0.09719623, 0.08526807,  0.07943536,  0.06797944,  0.09118502)
  tmbvcov[2,] <- c(0.04574391,  3.21994672,  0.91522073, 0.10980129, -0.28810783, -1.07845809,  0.51064309)
  tmbvcov[3,] <- c(0.09719623,  0.91522073,  0.40584816, 0.09981763, -0.01342937, -0.23826114,  0.21393310)
  tmbvcov[4,] <- c(0.08526807,  0.10980129,  0.09981763, 0.14821768,  0.05824110,  0.03110420,  0.08580658)
  tmbvcov[5,] <- c(0.07943536, -0.28810783, -0.01342937, 0.05824110,  0.16979550,  0.16423022,  0.02255625)
  tmbvcov[6,] <- c(0.06797944, -1.07845809, -0.23826114, 0.03110420,  0.16423022,  0.50464751, -0.10296956)
  tmbvcov[7,] <- c(0.09118502,  0.51064309,  0.21393310, 0.08580658,  0.02255625, -0.10296956,  0.22602059)
  expect_equal(summ$vcov[c(6:7, 1:5), c(6:7, 1:5)], tmbvcov, tol=1e-4)
  ## Stand error for sigma (original parameter)
  summ2 <- cmLaplace$summary(opt, originalScale = TRUE)
  expect_equal(summ2$params$stdError[2], 0.5472659, tol=1e-4)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() 
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
  ## TMB cpp code:
  #include <TMB.hpp>
  #template<class Type>
  #Type objective_function<Type>::operator() ()
  # {
  #   DATA_VECTOR(y);
  #   PARAMETER(mu);
  #   PARAMETER(sigmaTrans);
  #   PARAMETER_VECTOR(theta);
  #   // Transformation for sigma
  #   Type sigma = 100 * exp(sigmaTrans) / (1 + exp(sigmaTrans));
  #   // Negative log-likelihood
  #   Type ans = 0;
  #   vector<Type> p(5);
  #   for(int i = 0; i < 5; i++){
  #     p[i] = exp(theta[i]) / (1 + exp(theta[i]));
  #     ans -= dnorm(theta[i], mu, sigma, true) + dbinom(y[i], Type(10), p[i], true);
  #   }
  #   ADREPORT(sigma);
  #   return ans;
  # }
  ## TMB R code:
  # library(TMB)
  # compile("test.cpp")
  # dyn.load(dynlib("test"))
  # data <- list(y = m$y)
  # parameters <- list(mu = m$mu, sigmaTrans = logit(m$sigma/100), theta = m$theta)
  # ## Fit model
  # obj <- MakeADFun(data, parameters, random="theta", DLL="test")
  # tmbopt <- nlminb(obj$par, obj$fn, obj$gr)
  # tmbrep <- sdreport(obj, getJointPrecision = TRUE)
  # tmbvcov <- inverse(tmbrep$jointPrecision)
})

test_that("Laplace with no random effects (simple linear regression) works", {
  set.seed(1)
  x <- rnorm(5)
  y <- sapply(-1 + x, rnorm, n = 1, sd = 1)
  m <- nimbleModel(
    nimbleCode({
      a ~ dnorm(0, sd = 10.0)
      b ~ dnorm(0, sd = 10.0)
      sigma ~ dunif(0, 100)
      for(i in 1:5){
        mu_y[i] <- a + b*x[i]
        y[i] ~ dnorm(mu_y[i], sd = sigma)
      }
    }),
    constants = list(x = x),
    data = list(y = y),
    inits = list(a = -1, b = 1, sigma = 1),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  summ <- cmLaplace$summary(opt)
  ## Compare results with those from TMB
  expect_equal(opt$par, c(-0.8899436, 1.1940911, 0.5744841), tol = 1e-5)
  expect_equal(opt$value, -4.323288, tol = 1e-7)
  expect_equal(summ$params$stdError, c(0.2598061, 0.2988869, 0.1816661), tol = 1e-5)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() 
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt) # The warnings are normal
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
  
  ## TMB cpp code
  #include <TMB.hpp>
  #template<class Type>
  # Type objective_function<Type>::operator() ()
  # {
  #   DATA_VECTOR(y);
  #   DATA_VECTOR(x);
  #   PARAMETER(a);
  #   PARAMETER(b);
  #   PARAMETER(sigma);
  #   Type nll = -sum(dnorm(y, a+b*x, sigma, true));
  #   return nll;
  # }
  ## R code
  # compile("lm.cpp")
  # dyn.load(dynlib("lm"))
  # set.seed(1)
  # x <- rnorm(5)
  # y <- sapply(-1 + x, rnorm, n = 1, sd = 1)
  # data <- list(y=y, x=x)
  # parameters <- list(a=-1, b=1, sigma=1)
  # obj <- MakeADFun(data, parameters, DLL="lm")
  # obj$hessian <- TRUE
  # tmbres <- do.call("optim", obj)
  # tmbsumm <- summary(sdreport(obj))
})

test_that("Laplace with no priors for unconstrained parameters works", {
  ## Here we re-use some of tests above and remove priors for parameters
  ## Test 1
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd = 2)
      a ~ dnorm(mu, sd = 3)
      # mu ~ dnorm(0, sd = 5)
    }), data = list(y = 4), inits = list(a = -1),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, 4, tol = 1e-4) 
  expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))
  summ <- cmLaplace$summary(opt, originalScale = TRUE, calcRandomEffectsStdError = TRUE, returnJointCovariance = TRUE)
  expect_equal(summ$random$estimate, 4, tol = 1e-5)
  # Covariance matrix 
  vcov <- matrix(c(1/(1/4+1/9), 0, 0, 0), nrow = 2) + matrix(c(4/13, 1), ncol = 1) %*% (13) %*% t(matrix(c(4/13, 1), ncol = 1))
  expect_equal(vcov, summ$vcov, tol = 1e-6)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE()
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
  
  ## Test 2
  set.seed(1)
  x <- rnorm(5)
  y <- sapply(-1 + x, rnorm, n = 1, sd = 1)
  m <- nimbleModel(
    nimbleCode({
      sigma ~ dunif(0, 100)
      for(i in 1:5){
        mu_y[i] <- a + b*x[i]
        y[i] ~ dnorm(mu_y[i], sd = sigma)
      }
    }),
    constants = list(x = x),
    data = list(y = y),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))
  
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  summ <- cmLaplace$summary(opt)
  ## Compare results with those from TMB
  expect_equal(opt$par, c(0.5744841, -0.8899436, 1.1940911), tol = 1e-5)
  expect_equal(opt$value, -4.323288, tol = 1e-7)
  expect_equal(summ$params$stdError, c(0.1816661, 0.2598061, 0.2988869), tol = 1e-5)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() 
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt) 
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
  
  ## Test 3
  set.seed(1)
  y <- array(rnorm(8, 6, 5), dim = c(2, 2, 2)) 
  cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      # for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
      mu_a[1] <- 0.8 * mu[1]
      mu_a[2] <- 0.2 * mu[2]
      for(i in 1:2) a[i, 1:2] ~ dmnorm(mu_a[1:2], cov = cov_a[1:2, 1:2])
      for(i in 1:2) {
        for(j in 1:2) {
          y[1, j, i] ~ dnorm( 0.5 * a[i, 1], sd = 1.8) 
          y[2, j, i] ~ dnorm( 0.1 * a[i, 2], sd = 1.2)
        }
      }
    }),
    data = list(y = y),
    inits = list(a = matrix(c(-2, -3, 0,  -1), nrow = 2)),
    constants = list(cov_a = cov_a),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  
  expect_equal(opt$par, c(12.98392, 406.04878), tol = 1e-4)
  expect_equal(opt$value, -41.86976, tol = 1e-6)
  # Check covariance matrix
  summ <- cmLaplace$summary(opt, returnJointCovariance = TRUE)
  tmbvcov <- matrix(nrow = 6, ncol = 6)
  tmbvcov[1,] <- c(6.625000e+00, 4.687500e+00,  4.050000e+00,  4.050000e+00, -2.693817e-11, -2.695275e-11)
  tmbvcov[2,] <- c(4.687500e+00, 9.250000e+02,  2.965628e-11,  2.967848e-11,  1.800000e+02,  1.800000e+02)
  tmbvcov[3,] <- c(4.050000e+00, 2.951367e-11,  3.995242e+00,  2.484758e+00,  5.596302e-01, -5.596302e-01)
  tmbvcov[4,] <- c(4.050000e+00, 2.951367e-11,  2.484758e+00,  3.995242e+00, -5.596302e-01,  5.596302e-01)
  tmbvcov[5,] <- c(-2.691772e-11, 1.800000e+02,  5.596302e-01, -5.596302e-01,  3.684693e+01,  3.515307e+01)
  tmbvcov[6,] <- c(-2.691772e-11, 1.800000e+02, -5.596302e-01,  5.596302e-01,  3.515307e+01,  3.684693e+01)
  
  expect_equal(summ$vcov[c(5,6,1,3,2,4), c(5,6,1,3,2,4)], tmbvcov, tol = 1e-4)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() 
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
  
  ## Test 4
  m <- nimbleModel(
    nimbleCode({
      # for(i in 1:3) {
      #   mu[i] ~ dnorm(0, sd = 10)
      # }
      mu_a[1] <- mu[1] + mu[2]
      mu_a[2] <- mu[2] + mu[3]
      a[1] ~ dnorm(mu_a[1], sd = 2)
      y[1] ~ dnorm(a[1], sd = 3)
      a[2] ~ dnorm(mu_a[2], sd = 2)
      y[2] ~ dnorm(a[2], sd =3)
      y[3] ~ dnorm(mu[3], sd = 3)
    }),
    data = list(y = c(2, 3, 5)),
    inits = list(a = c(1, 2)),
    buildDerivs = TRUE
  )
  
  mLaplace <- buildLaplace(model = m, control = list(allowNonPriors = TRUE))
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE, allowNonPriors = TRUE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit
  
  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, c(4, -2, 5), tol = 1e-3)
  expect_equal(opt$value, -6.420377, tol = 1e-6)
  ## Check covariance matrix
  summ <- cmLaplace$summary(opt, returnJointCovariance = TRUE)
  
  ## Covariance matrix from TMB
  tmbvcov <- matrix(nrow = 5, ncol = 5)
  tmbvcov[1,] <- c( 35, -2.20000e+01,  9.000000e+00,  9.000000e+00, -9.000000e+00)
  tmbvcov[2,] <- c(-22,  2.20000e+01, -9.000000e+00,  8.463230e-13,  9.000000e+00)
  tmbvcov[3,] <- c( 9,  -9.00000e+00,  9.000000e+00, -3.462231e-13,  3.462231e-13)
  tmbvcov[4,] <- c( 9,   8.46323e-13, -3.462231e-13,  9.000000e+00,  3.462231e-13)
  tmbvcov[5,] <- c(-9,   9.00000e+00,  3.462231e-13,  3.462231e-13,  9.000000e+00)
  
  expect_equal(summ$vcov[c(3:5, 1:2), c(3:5, 1:2)], tmbvcov, tol=1e-5)
  
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() 
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)

})

test_that("Laplace with crossed random effects works", {
  library(lme4)
  data(Penicillin)
  N <- nrow(Penicillin) 
  plate <- rep(1:24, each = 6)
  np <- 24
  sample <- rep(1:6, 24)
  ns <- 6
  
  m <- nimbleModel(
    nimbleCode({
      ## Intercept
      beta ~ dnorm(0, sd = 100)
      ## Standard deviations
      sigma ~ dgamma(1.0, 1.0)
      sigma_p ~ dgamma(1.0, 1.0)
      sigma_s ~ dgamma(1.0, 1.0)
      ## Random effects for plate
      for(i in 1:np){
        mup[i] ~ dnorm(0, sd = sigma_p)
      }
      ## Random effects for sample
      for(i in 1:ns){
        mus[i] ~ dnorm(0, sd = sigma_s)
      }
      ## Observations
      for(i in 1:N){
        mu_y[i] <- beta + mus[sample[i]] + mup[plate[i]]
        y[i] ~ dnorm(mu_y[i], sd = sigma) 
      }
    }),
    constants = list(N = N, np = np, ns = ns, plate = plate, sample = sample),
    data = list(y = Penicillin$diameter),
    # inits = list(beta = 0, sigma = 1, sigma_p = 1, sigma_s = 1, mus = rep(0, ns), mup = rep(0, np)),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  opt <- cmLaplace$LaplaceMLE()
  nimres <- cmLaplace$summary(opt, calcRandomEffectsStdError = TRUE)
  
  lme4_fit <- lmer(diameter ~ 1 + (1|plate) + (1|sample), data = Penicillin, REML = FALSE)
  lme4res <- summary(lme4_fit)
  
  expect_equal(nimres$params$estimate[1], lme4res$coefficients[,"Estimate"], tol=1e-6)
  expect_equal(nimres$params$estimate[c(3,4,2)], as.data.frame(VarCorr(lme4_fit))[,"sdcor"], tol = 5e-5)
  expect_equal(nimres$params$stdError[1], lme4res$coefficients[,"Std. Error"], tol=5e-4)
  expect_equal(nimres$random$estimate[25:30], as.vector(t(ranef(lme4_fit)$sample)), tol = 5e-5)
  expect_equal(nimres$random$estimate[1:24], as.vector(t(ranef(lme4_fit)$plate)), tol = 5e-6)
})

test_that("Laplace with nested random effects works", {
  library(lme4)
  data(Pastes)
  lme4_fit <- lmer(strength ~ 1 + (1|batch) + (1|batch:cask), data = Pastes, REML = FALSE)
  lme4res <- summary(lme4_fit)
  
  m <- nimbleModel(
    nimbleCode({
      ## Intercept
      beta ~ dnorm(0, sd = 100)
      ## Standard deviations
      sigma ~ dgamma(1.0, 1.0)
      sigma1 ~ dgamma(1.0, 1.0)
      sigma2 ~ dgamma(1.0, 1.0)
      ## Random effects for batch
      for(i in 1:10){
        mub[i] ~ dnorm(0, sd = sigma1)
      }
      ## Random effects for batch:cask
      for(i in 1:30){
        mubc[i] ~ dnorm(0, sd = sigma2)
      }
      ## Observations
      for(i in 1:60){
        mu_y[i] <- beta + mub[batch[i]] + mubc[cask[i]]
        y[i] ~ dnorm(mu_y[i], sd = sigma) 
      }
    }),
    constants = list(batch = rep(1:10, each = 6), cask = rep(1:30, each = 2)),
    data = list(y = Pastes$strength),
    buildDerivs = TRUE
  )
  mLaplace <- buildLaplace(model = m)
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)
  ## It seems that default start values (0, 1, 1, 1) for this example do not work well 
  ## for optimisation; use c(2, 2, 2, 2) instead
  opt <- cmLaplace$LaplaceMLE(pStart = c(2,2,2,2))
  nimres <- cmLaplace$summary(opt, calcRandomEffectsStdError = TRUE)
  
  expect_equal(nimres$params$estimate[1], lme4res$coefficients[,"Estimate"], tol = 1e-5)
  expect_equal(nimres$params$estimate[c(4, 3, 2)], as.data.frame(VarCorr(lme4_fit))[,"sdcor"], tol = 5e-5)
  expect_equal(nimres$params$stdError[1], lme4res$coefficients[,"Std. Error"], tol = 5e-5)
  expect_equal(nimres$random$estimate[seq(1, 40, by = 4)], as.vector(t(ranef(lme4_fit)$batch)), tol = 5e-4)
  expect_equal(nimres$random$estimate[-seq(1, 40, by = 4)], as.vector(t(ranef(lme4_fit)$`batch:cask`)), tol = 5e-4)
})


# To do:
# Various structures of random effects groupings

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
