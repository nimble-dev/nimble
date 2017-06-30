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
    nimbleOptions(debugCppLineByLine = TRUE)
    nimbleOptions(pauseAfterWritingFiles = FALSE)
    nf <- nimbleFunction(
        name = 'example',
        run = function (arg1 = double(1)) 
        {
            out <- cosh(arg1)
            return(out)
            returnType(double(1))
        } 
    )
    cnf <- compileNimble(nf, dirName = file.path(Sys.getenv('HOME'), 'tmp'), projectName = 'tf',
                         control = list(debugCpp = TRUE))
    x <- matrix(1:4, 2, 2)
    expect_equal(nf(x), cnf(x))  ## FIXME cnf(x) is flattened to a vector.
})

test_that('Tensorflow backend works for basic math', {
    skip_if_not_installed('tensorflow')
    set.seed(0)
    sapply(testsVaried, test_math)
    sapply(testsBasicMath, test_math)
    sapply(testsReduction, test_math)
    sapply(testsComparison, test_math)
})

if (0)  ## These tests currently fail.
test_that('Tensorflow backend works for basic math', {
    skip_if_not_installed('tensorflow')
    set.seed(0)
    sapply(testsMoreMath, test_math)  ## Fails "probit/iprobit of vector".
    sapply(testsMatrix, test_math)  ## Fails "forwardsolve matrix-vector".
})
