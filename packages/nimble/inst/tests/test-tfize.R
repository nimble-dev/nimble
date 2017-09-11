source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))

context("Testing of Tensorflow-ization")

nimbleOptions(showCompilerOutput = TRUE)  # DO NOT SUBMIT

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
    nimbleOptions(experimentalEnableDerivs = TRUE)
    nimbleOptions(verboseErrors = TRUE)

    make_fun <- function() {
        nimbleFunction(
            setup = function(){},
            run = function(x = double(0)) {
                outList <- derivs(testMethod(x))
                returnType(ADNimbleList())
                return(outList)
            },
            methods = list(
                testMethod = function(x = double(0)) {
                    out <- log(x)
                    returnType(double(0))
                    return(out)
                }
            ),
            enableDerivs = list('testMethod')
        )()
    }
    fun1 <- make_fun()
    fun2 <- make_fun()
    temporarilyAssignInGlobalEnv(fun1)
    temporarilyAssignInGlobalEnv(fun2)

    # FIXME This fails with:
    # Error in nimType2TfDtype[[sym$type]] : 
    #   attempt to select less than one element in get1index
    # where the interim concatenate symbol is missing from symTab
    # (hence sym == NULL, sym$type == NULL).
    withTensorflowEnabled(
        tf_fun <- compileNimble(fun1, showCompilerOutput = TRUE, resetFunctions = TRUE)
    )
    c_fun <- compileNimble(fun2, showCompilerOutput = TRUE)
    
    Rderiv <- D(expression(log(x)), 'x')
    x <- 1.4
    expect_equal(c_fun$run(x)$gradient[1], eval(Rderiv))
    expect_equal(tf_fun$run(x)$gradient[1], eval(Rderiv))
})
