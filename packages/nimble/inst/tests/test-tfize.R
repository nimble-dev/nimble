source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))

context("Testing of Tensorflow-ization")

nimbleOptions(showCompilerOutput = TRUE)

test_that('Tensorflow implementation of axpy works', {
    temporarilyEnableTensorflow()

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

test_that('Tensorflow multiple statements', {
    temporarilyEnableTensorflow()
    nimbleOptions(debugCppLineByLine = TRUE)
    nimbleOptions(showCompilerOutput = TRUE)
    nimbleOptions(pauseAfterWritingFiles = FALSE)
    nimbleOptions(experimentalNewSizeProcessing = TRUE)

    nf <- nimbleFunction(
        name = 'example',
        run = function(arg1 = double(2))
        {
            ## This compiles into two separate tensorflow graphs, whereas
            ## we want it to be fused into a single graph.
            x <- arg1 %*% arg1
            out <- x %*% x
            return(out)
            returnType(double(2))
        } 
    )
    cnf <- compileNimble(nf, dirName = file.path(Sys.getenv('HOME'), 'tmp'), projectName = 'tf',
                         control = list(debugCpp = TRUE))
    x <- matrix(1:4, 2, 2)
    expect_equal(nf(x), cnf(x))
})

test_that('Tensorflow backend works for basic math', {
    temporarilyEnableTensorflow()

    set.seed(0)
    sapply(testsVaried, test_math, 'tensorflow')
    sapply(testsBasicMath, test_math, 'tensorflow')
    sapply(testsReduction, test_math, 'tensorflow')
    sapply(testsComparison, test_math, 'tensorflow')
    sapply(testsMoreMath, test_math, 'tensorflow')
    sapply(testsMatrix, test_math, 'tensorflow')
})
