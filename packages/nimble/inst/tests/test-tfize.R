source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))

context("Testing of Tensorflow-ization")

nimbleOptions(experimentalUseTensorflow = TRUE)
nimbleOptions(showCompilerOutput = TRUE)

test_that('Tensorflow implementation of axpy works', {
    skip_if_not_installed('tensorflow')
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

if (0)  ## Known failure.
test_that('Tensorflow example', {
    skip_if_not_installed('tensorflow')
    nf <- nimbleFunction(
        run = function(arg1 = double(2), arg2 = double(1)) 
        {
            out <- sd(arg1 %*% arg2)
            return(out)
            returnType(double(0))
        } 
    )
    cnf <- compileNimble(nf)
    x <- matrix(1:6, 2, 3)
    y <- 1:3
    expect_equal(nf(x, y), cnf(x, y))
})



if (0)  ## These math tests currently fail.
test_that('Tensorflow backend works for basic math', {
    skip_if_not_installed('tensorflow')
    set.seed(0)
    sapply(testsVaried, test_math)
    sapply(testsBasicMath, test_math)
})
