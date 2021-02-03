source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)

context("Testing of nimble options")

test_that('withNimbleOptions works for zero options', {
    expected <- 'foo'
    actual <- withNimbleOptions(list(), 'foo')
    expect_equal(expected, actual)
})

test_that('withNimbleOptions works for one option', {
    expected <- 'foo'
    nimbleOptions(verboseErrors = FALSE)
    expect_equal(getNimbleOption('verboseErrors'), FALSE)
    actual <- withNimbleOptions(list(verboseErrors = TRUE), {
        expect_equal(getNimbleOption('verboseErrors'), TRUE)
        'foo'
    })
    expect_equal(getNimbleOption('verboseErrors'), FALSE)
    expect_equal(expected, actual)
})

test_that('withNimbleOptions works for two options', {
    expected <- 'foo'
    nimbleOptions(verboseErrors = FALSE, compileOnly = FALSE)
    expect_equal(getNimbleOption('verboseErrors'), FALSE)
    expect_equal(getNimbleOption('compileOnly'), FALSE)
    actual <- withNimbleOptions(list(verboseErrors = TRUE, compileOnly = TRUE), {
        expect_equal(getNimbleOption('verboseErrors'), TRUE)
        expect_equal(getNimbleOption('compileOnly'), TRUE)
        'foo'
    })
    expect_equal(getNimbleOption('verboseErrors'), FALSE)
    expect_equal(getNimbleOption('compileOnly'), FALSE)
    expect_equal(expected, actual)
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
