source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of Tensorflow-ization")

nimbleOptions(experimentalUseTensorflow = TRUE)
nimbleOptions(showCompilerOutput = TRUE)

test_that('Tensorflow implementation of axpy works', {
    nf <- nimbleFunction(
        run = function(a = double(), x = double(1), y = double(1)) {
            z <- a * x + y
            return(z)
            returnType(double(1))
        })
    cnf <- compileNimble(nf)

    a <- 0.1
    x <- c(1, 2, 3)
    y <- c(4, 5, 6)
    expect_equal(nf(a, x, y), cnf(a, x, y))
})

## These math tests currently fail.
if (0) {
    source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))
    
    set.seed(0)
    ans1 <- sapply(testsVaried, test_math)    ## 12
    ans2 <- sapply(testsBasicMath, test_math) ## 70
}