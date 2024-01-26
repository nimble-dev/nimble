source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of nimIntegrate() in NIMBLE code")

test_that("Basic use of integrate in an RCfunction calling an RCfunction", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
      run = function(theta = double(0), lower = double(0), upper = double(0)) {
            # padding with a 0 is no longer necessary. But integrand should always take a vector.
            output <- integrate(integrand, lower, upper, theta)
            returnType(double(1))
            return(output)
        }
    )
    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)
    cfun <- compileNimble(fun)

    theta <- pi
    expected_result <- theta*0.375
    resultR <- fun(theta, .5, 1)
    resultC <- cfun(theta, .5, 1)
    expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
    expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
    expect_identical(resultC[3], 0, info = "unexpected error code")
})

test_that("Basic use of integrate with rel.tol and abs.tol expressions", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )

    fun <- nimbleFunction(
      run = function(theta = double(0), lower = double(0), upper = double(0)) {
        # padding with a 0 is no longer necessary. But integrand should always take a vector.
        abs.tol <- .00001
        rtvec <- c(0.000005, 0.000005)
        output <- integrate(integrand, lower, upper, theta, abs.tol = 2*abs.tol/(sqrt(2)^2), rel.tol = sqrt(sum(rtvec)^2))
            returnType(double(1))
            return(output)
        }
    )
    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)
    cfun <- compileNimble(fun)

    theta <- pi
    expected_result <- theta*0.375
    resultR <- fun(theta, .5, 1)
    resultC <- cfun(theta, .5, 1)
    expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
    expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
    expect_identical(resultC[3], 0, info = "unexpected error code")
})

test_that("Basic use of integrate in an RCfunction calling an RCfunction with vector args", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*sum(theta))
            returnType(double(1))
        }
    )

    fun <- nimbleFunction(
      run = function(theta = double(1), lower = double(0), upper = double(0)) {
            output <- integrate(integrand, lower, upper, theta)
            returnType(double(1))
            return(output)
        }
    )
    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)
    cfun <- compileNimble(fun)

    theta <- c(pi/2, pi/2)
    expected_result <- sum(theta)*0.375
    resultR <- fun(theta, .5, 1)
    resultC <- cfun(theta, .5, 1)
    expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
    expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
    expect_identical(resultC[3], 0, info = "unexpected error code")
})


test_that("Basic use of integrate in a nimbleFunction calling its own method", {
  nf <- nimbleFunction(
    setup=TRUE,
    run = function(theta = double(0), lower = double(0), upper = double(0)) {
                                        # padding with a 0 is no longer necessary. But integrand should always take a vector.
      output <- integrate(integrand, lower, upper, theta)
      returnType(double(1))
      return(output)
    },
    methods = list(
      integrand = function(x = double(1), theta = double(1)) {
        return(x*theta[1])
        returnType(double(1))
      }
    )
  )
  nf1 <- nf()

  cnf1 <- compileNimble(nf1)

  theta <- pi
  expected_result <- theta*0.375
  nf1$integrand # strange bug that integrate(integrand...) will not find integrand unless is has been invoked
  resultR <- nf1$run(theta, .5, 1)
  resultC <- cnf1$run(theta, .5, 1)
  expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
  expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
  expect_identical(resultC[3], 0, info = "unexpected error code")
})

test_that("Basic use of integrate in a nimbleFunction calling another nimbleFunction's method", {
  callee <- nimbleFunction(
    setup = TRUE,
    methods = list(
      integrand = function(x = double(1), theta = double(1)) {
        return(x*theta[1])
        returnType(double(1))
      }
    )
  )
  nf <- nimbleFunction(
    setup=function(c1) {},
    run = function(theta = double(0), lower = double(0), upper = double(0)) {
                                        # padding with a 0 is no longer necessary. But integrand should always take a vector.
      output <- integrate(c1$integrand, lower, upper, theta)
      returnType(double(1))
      return(output)
    }
  )
  c1 <- callee()
  nf1 <- nf(c1)

  cnf1 <- compileNimble(nf1)

  theta <- pi
  expected_result <- theta*0.375
  resultR <- nf1$run(theta, .5, 1)
  resultC <- cnf1$run(theta, .5, 1)
  expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
  expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
  expect_identical(resultC[3], 0, info = "unexpected error code")
})

# Following case is not supported (also not supported for nimOptim, I believe)
## test_that("Basic use of integrate in a nimbleFunction using a nimbleFunctionList", {
##   callee_base <- nimbleFunctionVirtual(
##     methods = list(integrand = function(x = double(1), theta = double(1)) {returnType(double(1))})
##   )
##   callee <- nimbleFunction(
##     contains = callee_base,
##     setup = TRUE,
##     methods = list(
##       integrand = function(x = double(1), theta = double(1)) {
##         return(x*theta[1])
##         returnType(double(1))
##       }
##     )
##   )
##   nf <- nimbleFunction(
##     setup=function(c1) {
##       callee_list <- nimbleFunctionList(callee_base)
##       callee_list[[1]] <- c1
##     },
##     run = function(theta = double(0), lower = double(0), upper = double(0)) {
##                                         # padding with a 0 is no longer necessary. But integrand should always take a vector.
##       output <- integrate(callee_list[[1]]$integrand, lower, upper, theta)
##       returnType(double(1))
##       return(output)
##     }
##   )
##   c1 <- callee()
##   nf1 <- nf(c1)

##   ncnf1 <- compileNimble(nf1)

##   theta <- pi
##   expected_result <- theta*0.375
##   resultR <- fun(theta, .5, 1)
##   resultC <- cfun(theta, .5, 1)
##   expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
##   expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
##   expect_identical(resultC[3], 0, info = "unexpected error code")
## })

test_that("basic use of integrate with literal values for bounds", {
    ## Using explicit values for bounds of `integrate`
    integrand2 <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun2 <- nimbleFunction(
        run = function(theta = double(0)) {
            output <- integrate(integrand2, .5, 1, theta)
            returnType(double(1))
            return(output)
        }
    )
    temporarilyAssignInGlobalEnv(integrand2)
    temporarilyAssignInGlobalEnv(fun2)
    cfun2 <- compileNimble(fun2)
    resultR <- fun2(theta)
    resultC <- cfun2(theta)
    expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
    expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
    expect_identical(resultC[3], 0, info = "unexpected error code")
})

test_that("Basic use in a nimbleFunction with infinity bound", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(dnorm(x))
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
        run = function(lower = double(0), upper = double(0)) {
            tmp <- rep(0,2)  # cannot be scalar
            output <- integrate(integrand, lower, upper, tmp)
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)
    cfun <- compileNimble(fun)

    resultR <- fun(0, Inf)
    resultC <- cfun(0, Inf)
    expect_equal(resultC[1], 0.5, info = "unexpected nimIntegrate result")
    expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
    expect_identical(resultC[3], 0, info = "unexpected error code")
})

test_that("Basic use of nimIntegrate with a bound that uses Eigen implementation", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(dnorm(x))
            returnType(double(1))
        }
    )

    fun <- nimbleFunction(
       run = function(lower = double(1), upper = double(0)) {
            output <- integrate(integrand, sum(lower*c(1,1)), upper, 0)
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)
    cfun <- compileNimble(fun)

    resultR <- fun(c(-1, 1), Inf)
    resultC <- cfun(c(-1, 1), Inf)
    expect_equal(resultC[1], 0.5, info = "unexpected nimIntegrate result")
    expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
    expect_identical(resultC[3], 0, info = "unexpected error code")
})

test_that("Error trapping", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(sin(x*theta[1]))
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
        run = function(theta = double(0), lower = double(0), upper = double(0)) {
            output <- integrate(integrand, lower, upper, theta)
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)
    cfun <- compileNimble(fun)

    theta <- 10000  # should be enough to prompt error
    expect_error(fun(theta, 0, 1), "maximum number of subdivisions reached")
    expect_error(cfun(theta, 0, 1), "integration error with code 1")


    integrand2 <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(sin(x*theta[1]))
            returnType(double(1))
        }
    )
    
    fun2 <- nimbleFunction(
        run = function(theta = double(0), lower = double(0), upper = double(0)) {
            tmp <- theta
            output <- integrate(integrand2, lower, upper, tmp, stop.on.error = FALSE)
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand2)
    temporarilyAssignInGlobalEnv(fun2)
    cfun2 <- compileNimble(fun2)

    theta <- 10000  # should be enough to prompt error
    resultR <- fun2(theta, 0, 1)
    resultC <- cfun2(theta, 0, 1)
    expect_identical(resultC[3], 1) 
})


test_that("Use in a user-defined model function", {
    code <- nimbleCode({
        mu ~ dnorm(0, sd = 100)
        y[1:n] ~ dintN(mu)
    })

    integrand <- nimbleFunction(
        run = function(sigma = double(1), pars = double(1)) {
                                        # pars is (mu,y1,...,yn)
            n <- length(pars)-1
            y <- pars[2:(n+1)]
            ybar <- mean(y)
            ss <- sum((y-ybar)^2)
            result <- -(n/2)*log(2*pi) -n*log(sigma) -0.5 * (ss + n*(ybar-pars[1])^2)/sigma^2
            return(exp(result))
            returnType(double(1))
        }
    )

    dintN <- nimbleFunction(
        run = function(x = double(1), mu = double(0),
                       log = integer(0, default = 0)) {
            returnType(double(0))
            pars <- c(mu,x)
            ## 0,5 imply uniform prior on sigma on (0,5)
            prob <- integrate(integrand, 0, 5, pars)[1]
            if(log) return(log(prob)) else return(prob)
        })

    rintN <- nimbleFunction(
        run=function(n = integer(0), mu = double(0)) {
            returnType(double(1))
            return(rep(0, 2))
        })
    
    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(dintN)
    temporarilyAssignInGlobalEnv(rintN)

    set.seed(1)
    y  <- rnorm(5)
    m <- nimbleModel(code, data = list(y = y), constants = list(n=5),
                     inits = list(mu = 0))
    cm <- compileNimble(m)
    expect_silent(value <- cm$calculate('y'))

    mcmc <- buildMCMC(m)
    cmcmc <- compileNimble(mcmc,project=m)
    expect_silent(out <- runMCMC(cmcmc, 50))
   
})

test_that("Error trapping for invalid types for `params` or non-scalar `lower`, `upper`", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
        run = function(lower = double(0), upper = double(0)) {
            output <- integrate(integrand, lower, upper)
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)
    expect_error(cfun <- compileNimble(fun), "argument must be provided")

    integrand2 <- nimbleFunction(
        run = function(x = double(1), theta = double(0)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun2 <- nimbleFunction(
        run = function(theta = double(2), lower = double(0), upper = double(0)) {
            output <- integrate(integrand2, lower, upper, theta)
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand2)
    temporarilyAssignInGlobalEnv(fun2)
    expect_error(cfun2 <- compileNimble(fun2), "must be a one-dimensional array")   

    integrand3 <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun3 <- nimbleFunction(
        run = function(theta = double(2), lower = double(0), upper = double(0)) {
            output <- integrate(integrand3, lower, upper, theta + 3)
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand3)
    temporarilyAssignInGlobalEnv(fun3)
    expect_error(cfun3 <- compileNimble(fun3), "must be a one-dimensional array")   

    integrand4 <- nimbleFunction(
        run = function(x = double(1), theta = double(0)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun4 <- nimbleFunction(
        run = function(theta = double(0), lower = double(0), upper = double(0)) {
            output <- integrate(integrand4, c(lower,1), upper, c(theta, 0))
            returnType(double(1))
            return(output)
        }
    )

    temporarilyAssignInGlobalEnv(integrand4)
    temporarilyAssignInGlobalEnv(fun4)
    expect_error(cfun4 <- compileNimble(fun4), "must be scalars")   


    integrand5 <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun5 <- nimbleFunction(
        run = function(theta = double(0), lower = double(0), upper = double(0)) {
            tmp <- c(lower, 0)
            output <- integrate(integrand5, tmp, upper, 0)
            returnType(double(1))
            return(output)
        }
    )
    temporarilyAssignInGlobalEnv(integrand5)
    temporarilyAssignInGlobalEnv(fun5)
    expect_error(cfun5 <- compileNimble(fun5), "must be scalars")   
})

test_that("integrate with `param` expressions works", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1]/theta[2])
            returnType(double(1))
        }
    )

    fun <- nimbleFunction(
        run = function(theta = double(1), lower = double(0), upper = double(0)) {
            output <- integrate(integrand, lower, upper, exp(theta))
            returnType(double(1))
            return(output)
        }
    )
    temporarilyAssignInGlobalEnv(integrand)
    temporarilyAssignInGlobalEnv(fun)

    cfun <- compileNimble(fun)
    theta <- c(2*pi, 2)
    expected_result <- exp(theta[1])*0.375/exp(theta[2])
    resultR <- fun(theta, .5, 1)
    resultC <- cfun(theta, .5, 1)
    expect_equal(resultC[1], expected_result, info = "unexpected nimIntegrate result")
    expect_equal(resultR[1], resultC[1], info = "disparity in compiled and uncompiled nimIntegrate result")
    expect_identical(resultC[3], 0, info = "unexpected error code")
})
