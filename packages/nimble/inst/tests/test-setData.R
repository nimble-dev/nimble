source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

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

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
