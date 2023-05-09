source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

# The methods "SANN" and "Brent" are not supported.
methodsAllowingGradient <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")
methodsAllowingBounds <- c("L-BFGS-B")

# See test-optim for non-AD optim tests

test_that("optim() respects parscale in R and C++ while using AD", {
  nf <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      return(sum(x^2))
      returnType(double())
    },
    methods = list(
      gradRun = function(x = double(1)) {
        d <- derivs(run(x), wrt = 1:length(x), order = 1)
        ans <- d$jacobian[1,]
        return(ans)
        returnType(double(1))
      },
      optRun = function(x = double(1), parscale = double(1),
                        method = character()) {
        con <- nimOptimDefaultControl()
        con$parscale <- parscale
        ans <- optim(x, run, gr = gradRun, method = method, control = con)
        return(ans)
        returnType(optimResultNimbleList())
      }
    ),
    buildDerivs = 'run'
  )
  nf1 <- nf()
  cnf1 <- compileNimble(nf1)
  junk <- nf1$run # needed to bring run into scope (wierd)
  junk <- nf1$gradRun # ditto
  for(method in c(methodsAllowingGradient, methodsAllowingBounds)) {
    optR <- nf1$optRun(c(1, 1), c(1, 1), method = method)
    optC <- cnf1$optRun(c(1, 1), c(1, 1), method = method)
    expect_equal(optR$par, optC$par)

    optR <- nf1$optRun(c(1, 1), c(3,2), method = method)
    optC <- cnf1$optRun(c(1, 1), c(3,2), method = method)
    expect_equal(optR$par, optC$par)
  }
})

test_that("optim() respects parscale for Hessian in R and C++ using AD gradients", {
  nf <- nimbleFunction(
    setup = TRUE,
    run = function(x = double(1)) {
      xp <- c(x[1] / 2, x[2] / 4)
      return(sum(xp^2) + 0.5*prod(xp))
      returnType(double())
    },
    methods = list(
      gradRun = function(x = double(1)) {
        d <- derivs(run(x), wrt = 1:length(x), order = 1)
        ans <- d$jacobian[1,]
        return(ans)
        returnType(double(1))
      },
      optRun = function(x = double(1), parscale = double(1),
                        method = character()) {
        con <- nimOptimDefaultControl()
        con$parscale <- parscale
        ans <- optim(x, run, gr = gradRun, method = method, control = con, hessian=TRUE)
        return(ans)
        returnType(optimResultNimbleList())
      }
    ),
    buildDerivs = 'run'
  )
  nf1 <- nf()
  cnf1 <- compileNimble(nf1)
  junk <- nf1$run # needed to bring run into scope (wierd)
  junk <- nf1$gradRun # ditto
  for(method in c(methodsAllowingGradient, methodsAllowingBounds)) {
    optR <- nf1$optRun(c(1, 1), c(1, 1), method = method)
    optC <- cnf1$optRun(c(1, 1), c(1, 1), method = method)
    expect_equal(optR$par, optC$par)
    expect_equal(optR$hessian, optC$hessian)

    optR <- nf1$optRun(c(1, 1), c(3,2), method = method)
    optC <- cnf1$optRun(c(1, 1), c(3,2), method = method)
    expect_equal(optR$par, optC$par)
    expect_equal(optR$hessian, optC$hessian)
  }
})
