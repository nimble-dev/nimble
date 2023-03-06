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
                                              methods = 1:3 # methods to check
                                              ) {
  vars <- cm$getVarNames()
  reset <- function() {
    for(v in vars) cm[[v]] <- m[[v]]
  }
  if(missing(opt)) {
    reset()
    opt <- cL$LaplaceMLE()
  }
  ref_method <- cL$get_method()
  for(method in methods) {
    if(method != ref_method) {
      reset()
      cL$set_method(method)
      opt_alt <- cL$LaplaceMLE()
      expect_equal(opt$par, opt_alt$par, tolerance = 0.01)
      expect_equal(opt$value, opt_alt$value, tolerance = 1e-7)
    }
  }
  invisible(NULL)
}

test_that("Laplace simplest 1D works", {
  m <- nimbleModel(
    nimbleCode({
      y ~ dnorm(a, sd=2)
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

  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, 4, tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y] = 9 + 4 = 13
  # Cov[a, y] = V[a] = 9 (not needed)
  expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))

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

  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-2)

  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)

  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
})

# Following currently fails because getConditionallyIndependentSets does not trace upward to stochastic nodes
test_that("Laplace with 3x2D random effects for 1D data needing joint integration works, with intermediate nodes", {
  set.seed(1)
  # y[i, j] is jth datum from ith group
  y <- array(rnorm(8, 6, 5), dim = c(2, 2, 2)) #

  cov_a <- matrix(c(2, 1.5, 1.5, 2), nrow = 2)
  m <- nimbleModel(
    nimbleCode({
      for(i in 1:2) mu[i] ~ dnorm(0, sd = 10)
      mu_a[1] <- 0.8 * mu[1]
      mu_a[2] <- 0.2 * mu[2]
      for(i in 1:2) a[i, 1:2] ~ dmnorm( mu_a[1:2], cov = cov_a[1:2, 1:2])
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

  # For this case, we build up the correct answer more formulaically
  # Define A as the vector a[1, 1], a[1, 2], a[2, 1], a[2, 2]
  cov_A <- matrix(0, nrow = 4, ncol = 4)
  cov_A[1:2, 1:2] <- cov_a
  cov_A[3:4, 3:4] <- cov_a
  # Define Y as the vector y[1,1,1],y[2,1,1],y[1,2,1],y[2,2,1], then same with last index 2
  # Define E[Y] as IA %*% A, where:
  IA <- matrix(0, nrow = 8, ncol = 4)
  IA[c(1, 3), 1] <- 0.5
  IA[c(2, 4), 2] <- 0.1
  IA[c(5, 7), 3] <- 0.5
  IA[c(6, 8), 4] <- 0.1

  # define cov_y_given_a as the Cov[Y | A]
  cov_y_given_a <- matrix(0, nrow = 8, ncol = 8)
  diag(cov_y_given_a) <- rep(c(1.8^2, 1.2^2), 4)
  # And finally get cov_Y, the marginal (over A) covariance of Y
  cov_Y <- IA %*% cov_A %*% t(IA) + cov_y_given_a
  chol_cov <- chol(cov_Y)

  # make a log likelihood function
  nlogL <- function(mu) {
    mean_Y <- rep(c(0.8*0.5*mu[1], 0.2*0.1*mu[2]), 4)
    -dmnorm_chol(as.numeric(y), mean_Y, cholesky = chol_cov, prec_param=FALSE, log = TRUE)
  }
  # maximize it
  opt_manual <- optim(c(20, 100), nlogL, method = "BFGS")

  expect_equal(opt$par, opt_manual$par, tol = 1e-4)
  expect_equal(opt$value, -opt_manual$value, tol = 1e-5)
  for(v in cm$getVarNames()) cm[[v]] <- m[[v]]
  optNoSplit <- cmLaplaceNoSplit$LaplaceMLE() # some warnings are ok here
  expect_equal(opt$par, optNoSplit$par, tol = 1e-4)
  expect_equal(opt$value, optNoSplit$value, tol = 1e-7)
  check_laplace_alternative_methods(cmLaplace, cm, m, opt)
  check_laplace_alternative_methods(cmLaplaceNoSplit, cm, m, optNoSplit)
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
    constants = list(g = g, ng = max(g), n = n, x = x)
  )
  params <- c("fixed_int", "fixed_slope", "sigma_int", "sigma_slope", "sigma_res")
  values(m, params) <- c(10, 0.5, 3, .25, 0.2)
  m$simulate(m$getDependencies(params, self = FALSE))
  m$setData('y')
  y <- m$y
  library(lme4)

  manual_fit <- lmer(y ~ x + (1 + x || g), REML = FALSE)

  mLaplace <- buildLaplace(model = m)
  mLaplaceNoSplit <- buildLaplace(model = m, control = list(split = FALSE))
  cm <- compileNimble(m)
  cL <- compileNimble(mLaplace, mLaplaceNoSplit, project = m)
  cmLaplace <- cL$mLaplace
  cmLaplaceNoSplit <- cL$mLaplaceNoSplit

  opt <- cmLaplace$LaplaceMLE()
  # This test works. The opt$par matches the manual_fit.
  # But it is still obscured by the parameter transformation,
  # and while allowing that as an option, we also need to
  # work how to handle a summaryLaplace feature.
})



# To do:
# Multivariate random effects that are separable
# Multiple parameters
# parameters needing transformation
# Various structures of random effects groupings
# Crossed scalar random effects
# Correlated intercept and slope
# Try having no prior

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
