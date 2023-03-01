# Tests of Laplace approximation
source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'testthat', 'AD_test_utils.R'), package = 'nimble'))
EDopt <- nimbleOptions("enableDerivs")
BMDopt <- nimbleOptions("buildModelDerivs")
nimbleOptions(enableDerivs = TRUE)
nimbleOptions(buildModelDerivs = TRUE)
nimbleOptions(allowDynamicIndexing = FALSE)

library(nimble)

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
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)

  opt <- cmLaplace$LaplaceMLE()
  expect_equal(opt$par, 4, tol = 1e-4) # optim's reltol is about 1e-8 but that is for the value, not param.
  # V[a] = 9
  # V[y] = 9 + 4 = 13
  # Cov[a, y] = V[a] = 9 (not needed)
  expect_equal(opt$value, dnorm(4, 4, sd = sqrt(13), log = TRUE))
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
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)

  opt <- cmLaplace$LaplaceMLE() # some warnings are ok here.
  expect_equal(opt$par, 40, tol = 1e-4) # 40 = 4 * (1/.2) * (1/.5)
  # V[a] = 9
  # V[y] = 0.2^2 * 9 + 4 = 4.36
  expect_equal(opt$value, dnorm(0.1*40, 0.1*40, sd = sqrt(4.36), log = TRUE))
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
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)

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
  cm <- compileNimble(m)
  cmLaplace <- compileNimble(mLaplace, project = m)

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
})

# To do:
# Multivariate random effects
# Multiple parameters
# Various structures of random effects groupings

nimbleOptions(enableDerivs = EDopt)
nimbleOptions(buildModelDerivs = BMDopt)
