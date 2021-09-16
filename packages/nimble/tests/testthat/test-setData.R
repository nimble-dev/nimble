source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

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

nimbleOptions(verbose = TRUE)

test_that("ignores extra variable not in model", {
    expect_message(model$setData('a', 'b', 'c'),
                   'is not a variable in the model')
    model$resetData()
})

test_that("unnamed argument produces error", {
    expect_warning(model$setData(a = c(11, NA, 13:15), 7), 'unnamed element')
    model$resetData()
})

nimbleOptions(verbose = FALSE)

test_that("unnamed arguments produce error", {
    expect_error(model$setData('a', 3), 'multiple inputs must be named')
})

test_that("no named arguments produces error", {
    expect_error(model$setData(3, 7), 'multiple inputs must be named')
})

test_that("weird character input produces error", {
    expect_error(model$setData(c('a','d'),'b'))
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

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
