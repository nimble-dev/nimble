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

if (0)  # Known failure.
test_that('Tensorflow implementation can compile parentheses', {
    nf <- nimbleFunction(
        run = function(arg1 = double(2), arg2 = double(1)) 
        {
            out <- sd(arg1 %*% arg2)
            return(out)
            returnType(double(0))
        } 
    )
    cnf <- compileNimble(nf)
    expect_equal(nf(x), cnf(x))
})


## These math tests currently fail.
if (0) {
    source(system.file(file.path('tests', 'mathTestLists.R'), package = 'nimble'))
    
    set.seed(0)
    ans1 <- sapply(testsVaried, test_math)    ## 12
    ans2 <- sapply(testsBasicMath, test_math) ## 70
}
