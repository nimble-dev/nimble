source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))

context("Testing of Tensorflow-ization")

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

test_that('Tensorflow works with nimDerivs()', {
    temporarilyEnableTensorflow()
    nimbleOptions(debugCppLineByLine = TRUE)
    nimbleOptions(showCompilerOutput = TRUE)
    nimbleOptions(pauseAfterWritingFiles = FALSE)
    nimbleOptions(experimentalNewSizeProcessing = TRUE)
    
    nf1 <- nimbleFunction(
        name = 'example',
        run = function(x = double(1)) {
            return(sum(x))
            returnType(double(0))
        }
    )
    nf2 <- nimbleFunction(
        name = 'example_gradient',
        run = function(arg1 = double(1)) {
            return(nimDerivs(nf1, order = c(0, 1)))
            returnType(ADNimbleList())
        }
    )
    dirName <- file.path(Sys.getenv('HOME'), 'tmp')
    cnf1 <- compileNimble(nf1, dirName = dirName, projectName = 'tf',
                          control = list(debugCpp = TRUE))
    cnf2 <- compileNimble(nf2, dirName = dirName, projectName = 'tf',
                          control = list(debugCpp = TRUE))
    x <- c(-2, -1, 0, 1, 2)
    dx <- c(-4, -2, 0, 2, 4)
    expect_equal(cnf(x), dx)
})

test_that('Tensorflow can compute derivatives', {
    nimbleOptions(experimentalEnableDerivs = TRUE,
                  experimentalUseTensorflow = FALSE)

    ADfun <- nimbleFunction(
        setup = function(){},
        run = function(x = double(0)) {
            outList <- derivs(testMethod(x))
            returnType(ADNimbleList())
            return(outList)
        },
        methods = list(
            testMethod = function(x = double(0)) {
                out <- dnorm(x,0,1)
                returnType(double())
                return(out)
            }
        ),
        enableDerivs = list('testMethod')
    )
    
    ADfunInst <- ADfun()
    temporarilyAssignInGlobalEnv(ADfunInst)
    cADfunInst <- compileNimble(ADfunInst, showCompilerOutput = TRUE)
    
    Rderiv <- D(expression((1/(sqrt(2*pi)))*exp(-(x^2)/2)), 'x') ## Can be replaced by NIMBLE's R version of derivs in the future.
    x <- 1.4
    expect_equal(cADfunInst$run(x)$gradient[1], eval(Rderiv)) ## Temporary simple test to make sure compilation and gradient calculation work.
})
