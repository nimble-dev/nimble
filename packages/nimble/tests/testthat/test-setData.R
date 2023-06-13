source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = TRUE)

context("Testing setData")

model <- nimbleModel(
    nimbleCode({
        for(i in 1:5) {
            a[i] ~ dnorm(0,1)
            b[i] ~ dnorm(0,1)
        }
    }))

test_that("set by current value and try NA", {
    model$a <- c(1, NA, 3:5)
    model$setData('a')
    expect_identical(as.numeric(model$a), as.numeric(c(1, NA, 3:5)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
})

test_that("set by new value and try NA", {
    model$a <- 1:5
    model$setData(a = c(6, NA, 8:10))
    expect_identical(as.numeric(model$a), as.numeric(c(6, NA, 8:10)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
})

test_that("set by new value from list and try NA", {
    model$a <- 1:5
    model$setData(list(a = c(6, NA, 8:10)))
    expect_identical(as.numeric(model$a), as.numeric(c(6, NA, 8:10)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
})

test_that("set by new value from data frame and try NA", {
    model$a <- 1:5
    model$setData(data.frame(a = c(6, NA, 8:10)))
    expect_identical(as.numeric(model$a), as.numeric(c(6, NA, 8:10)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
})

test_that("set both to current values by character vector", {
    model$a <- c(1, NA, 3:5)
    model$b <- c(6, NA, 8:10)
    model$setData(c('a','b'))
    expect_identical(as.numeric(model$a), as.numeric(c(1, NA, 3:5)))
    expect_identical(as.numeric(model$b), as.numeric(c(6, NA, 8:10)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    expect_identical(model$isData('b'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
    expect_equal(all(model$isData('b')), FALSE)
})

test_that("set both to current values by two character args", {
    model$a <- c(1, NA, 3:5)
    model$b <- c(6, NA, 8:10)
    model$setData('a','b')
    expect_identical(as.numeric(model$a), as.numeric(c(1, NA, 3:5)))
    expect_identical(as.numeric(model$b), as.numeric(c(6, NA, 8:10)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    expect_identical(model$isData('b'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
    expect_equal(all(model$isData('b')), FALSE)
})

test_that("set one to current values and one to new value by list", {
    model$a <- c(1, NA, 3:5)
    model$b <- c(6, NA, 8:10)
    model$setData(list(a = c(11, NA, 13:15),'b'))
    expect_identical(as.numeric(model$a), as.numeric(c(11, NA, 13:15)))
    expect_identical(as.numeric(model$b), as.numeric(c(6, NA, 8:10)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    expect_identical(model$isData('b'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
    expect_equal(all(model$isData('b')), FALSE)
})

test_that("set one to current values and one to new value by two arguments", {
    model$a <- c(1, NA, 3:5)
    model$b <- c(6, NA, 8:10)
    model$setData(a = c(11, NA, 13:15),'b')
    expect_identical(as.numeric(model$a), as.numeric(c(11, NA, 13:15)))
    expect_identical(as.numeric(model$b), as.numeric(c(6, NA, 8:10)))
    expect_identical(model$isData('a'), c(TRUE, FALSE, rep(TRUE, 3)))
    expect_identical(model$isData('b'), c(TRUE, FALSE, rep(TRUE, 3)))
    model$resetData()
    expect_equal(all(model$isData('a')), FALSE)
    expect_equal(all(model$isData('b')), FALSE)
})

test_that("ignores extra variable not in model", {
    expect_message(model$setData('a', 'b', 'c'),
                   'is not a variable in the model')
    model$resetData()
})

test_that("unnamed argument produces error", {
    expect_warning(model$setData(a = c(11, NA, 13:15), 7), 'unnamed element')
    model$resetData()
})

test_that("unnamed arguments produce error", {
    expect_error(model$setData('a', 3), 'multiple inputs must be named')
})

test_that("no named arguments produces error", {
    expect_error(model$setData(3, 7), 'multiple inputs must be named')
})

test_that("weird character input produces error", {
    expect_error(model$setData(c('a','d'),'b'))
})

test_that("deterministic data nodes produce warning", {
    expect_message(
        model <- nimbleModel(
            nimbleCode({
                for(i in 1:5) {
                    z[i] <- mu[i]
                }
                z[6] ~ dnorm(0, 1)
                y ~ dnorm(0, 1)
            }), data = list(z = rnorm(6), y = rnorm(1))),
        "Deterministic node values will not be flagged as 'data'"
    )
    expect_identical(model$isData('z', 1), c(rep(FALSE, 5), TRUE))

    expect_message(
        model <- nimbleModel(
            nimbleCode({
                for(i in 1:5) {
                    z[i] <- mu[i]
                }
                z[6] ~ dnorm(0, 1)
                y ~ dnorm(0, 1)
            }), constants = list(z = rnorm(6))),
        "Deterministic node values will not be flagged as 'data'"
    )

    expect_message(
        model <- nimbleModel(
            nimbleCode({
                z[1:2] <- mu[1:2]
                z[3] ~ dnorm(0, 1)
                y ~ dnorm(0, 1)
            }), data = list(z = rnorm(3), y = rnorm(1))),
        "Deterministic node values will not be flagged as 'data'"
    )

    expect_message(
        model <- nimbleModel(
            nimbleCode({
                z[1:2] <- mu[1:2]
                z[3] ~ dnorm(0, 1)
                y ~ dnorm(0, 1)
            }), data = list(z = c(rep(NA, 2), rnorm(1)), y = rnorm(1))),
        "Checking model sizes"
    )
    
    expect_message(
        model <- nimbleModel(
            nimbleCode({
                for(i in 1:5) {
                    z[i] <- mu[i]
                }
                z[6] ~ dnorm(0, 1)
                y ~ dnorm(0, 1)
            }), data = list(z = c(rep(NA, 5), rnorm(1)), y = rnorm(1))),
        "Checking model sizes"        
    )
})

test_that("RHSonly nodes not flagged as data", {
    model <- nimbleModel(
        nimbleCode({
            y ~ dnorm(mu, 1)
        }), data = list(y = 1, mu = 2))
    expect_false(model$isData('mu'))

    model <- nimbleModel(
        nimbleCode({
            y ~ dnorm(mu, 1)
        }), data = list(y = 1), constants = list(mu = 2))
    expect_length(model$isData('mu'), 0)
})

test_that("data frame as matrix processed correctly", {
    model <- nimbleModel(
        nimbleCode({
            for(i in 1:5) 
                for(j in 1:2) 
                    a[i,j] ~ dnorm(0,1)
        }))
    v1vals <- rnorm(5)
    v2vals <- rnorm(5)
    model$setData(a = data.frame(v1 = v1vals, v2 = v2vals))
    tmp <- as.matrix(data.frame(v1 = v1vals, v2 = v2vals))
    dimnames(tmp) <- NULL
    expect_identical(model$a, tmp)
    expect_equal(all(model$isData('a')), TRUE)
    model$resetData()
    expect_error(model$setData(a = data.frame(v1 = c('a','b','c','d','e'), v2 = v2vals)),
                 'must be numeric')
})

model <- nimbleModel(
    nimbleCode({
        for(i in 1:1)
            a[i] ~ dnorm(0,1)
    }))

test_that("set 1-length data vector", {
    expect_silent(model$setData(a = 3))
    expect_equal(model$isData('a'), TRUE)
})

test_that("mixed data/NAs in multivariate nodes treated as a data node", {
    code <- nimbleCode({
        y[1:3] ~ dmnorm(mu[1:3], pr[1:3,1:3])
        y[4] ~ dnorm(0,1)
    })
    m  <- nimbleModel(code)
    ## first element missing
    m$setData(y = c(NA, 3, 5, 3))
    expect_true(m$isData('y[4]'))
    expect_true(m$isData('y[1:3]'))
    expect_identical(m$isDataEnv$y, c(FALSE, TRUE, TRUE, TRUE))
    ## non-first element missing
    m$setData(y = c(3, NA, 5, NA))
    expect_false(m$isData('y[4]'))
    expect_true(m$isData('y[1:3]'))
    expect_identical(m$isDataEnv$y, c(TRUE, FALSE, TRUE, FALSE))
})

test_that("mixed data and non-data in variable with 'missing' nodes" {
    ## This is meant to mimic capture-recapture type situations,
    ## and test bug in issue 1326 (and issue 1324).
    set.seed(1)
    code <- nimbleCode({
        for(i in 1:3)
            for(j in (start[i]+1):end[i])
                mu[i,j] ~ dnorm(mu[i,j-1],1)
        mu[1,1] <- 7  # add deterministic to check that too
    })
    start <- c(3,1,5)
    end <- c(5,2,8)

    data <- list(mu = matrix(rnorm(3*8),3))
    data$mu[1,5] <- NA
    data$mu[2,1] <- NA
    inits <- list(mu = matrix(rnorm(3*8),3))
    expect_message(m <- nimbleModel(code, constants = list(start = start, end = end),
                     data = data, inits = inits, calculate = FALSE), "Ignoring non-NA values in `inits`")

    expected <- matrix(TRUE, 3, 8)
    ## determ and RHSonly are not data, nor is anyting initialized with NA in 'data'.
    expected[1,3] <- expected[2,1] <- expected[3,5] <- expected[1,1] <- expected[1,5] <- FALSE
    expect_identical(m$isDataEnv$mu, expected)
    expected <- data$mu
    ## Only NA values in 'data' will be overwritten.
    expected[1,5] <- inits$mu[1,5]
    expected[2,1] <- inits$mu[2,1]
    expect_identical(m$mu, expected)

})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
