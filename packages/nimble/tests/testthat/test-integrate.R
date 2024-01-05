source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of nimIntegrate() in NIMBLE code")

test_that("Basic use in a nimbleFunction", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
        run = function(theta = double(0), lower = double(0), upper = double(0)) {
            tmp = c(theta, 0)  # cannot be scalar, so pad with zero.
            output = integrate(integrand, lower, upper, tmp)
            returnType(double(1))
            return(output)
        }
    )

    cfun <- compileNimble(fun)

    theta <- pi
    expected_result <- theta*0.375
    resultR <- fun(theta, .5, 1)
    resultC <- cfun(theta, .5, 1)
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
            tmp = rep(0,2)  # cannot be scalar
            output = integrate(integrand, lower, upper, tmp)
            returnType(double(1))
            return(output)
        }
    )

    cfun <- compileNimble(fun)

    resultR <- fun(0, Inf)
    resultC <- cfun(0, Inf)
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
            tmp = c(theta, 0)  # cannot be scalar, so pad with zero.
            output = integrate(integrand, lower, upper, tmp)
            returnType(double(1))
            return(output)
        }
    )

    cfun <- compileNimble(fun)

    theta <- 10000  # should be enough to prompt error
    expect_error(fun(theta, 0, 1), "maximum number of subdivisions reached")
    expect_error(cfun(theta, 0, 1), "integration error with code 1")


    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(sin(x*theta[1]))
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
        run = function(theta = double(0), lower = double(0), upper = double(0)) {
            tmp = c(theta, 0)  # cannot be scalar, so pad with zero.
            output = integrate(integrand, lower, upper, tmp, stop.on.error = FALSE)
            returnType(double(1))
            return(output)
        }
    )

    cfun <- compileNimble(fun)

    theta <- 10000  # should be enough to prompt error
    resultR <- fun(theta, 0, 1)
    resultC <- cfun(theta, 0, 1)
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
                                        # 0,5 imply uniform prior on sigma on (0,5)
            prob <- integrate(integrand, 0, 5, pars)
            if(log) return(log(prob)) else return(prob)
        })

    set.seed(1)
    y  <- rnorm(5)
    m <- nimbleModel(code, data = list(y = y), constants = list(n=5),
                     inits = list(mu = 0))
    cm <- compileNimble(m)
    cm$calculate('y')

    mcmc <- buildMCMC(m)
    cmcmc <- compileNimble(mcmc,project=m)
    out <- runMCMC(cmcmc, 5000))
   
}



test_that("Error trapping for non-vector 'params'", {
    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(1)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
        run = function(lower = double(0), upper = double(0)) {
            output = integrate(integrand, lower, upper)
            returnType(double(1))
            return(output)
        }
    )

    expect_error(cfun <- compileNimble(fun), "argument must be provided")

    integrand <- nimbleFunction(
        run = function(x = double(1), theta = double(0)) {
            return(x*theta[1])
            returnType(double(1))
        }
    )
    
    fun <- nimbleFunction(
        run = function(theta = double(0), lower = double(0), upper = double(0)) {
            output = integrate(integrand, lower, upper, theta)
            returnType(double(1))
            return(output)
        }
    )

    expect_error(cfun <- compileNimble(fun), "must be a one-dimensional array")   
})
