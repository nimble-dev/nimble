source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of the optim() function in NIMBLE code")

test_that("nimbleFunction of optim stub works", {
    nimFun <- nimbleFunction(
        run = function(x = double(0)) {
            return(optim(x))
            returnType(double(0))
        }
    )
    ans <- nimFun(0)
    expect_equal(ans, 0)
})

test_that("compileNimble of optim stub works", {
    nimFun <- nimbleFunction(
        run = function(x = double(0)) {
            return(optim(x))
            returnType(double(0))
        }
    )
    compiledFun <- compileNimble(nimFun)
    ans <- compiledFun$run(0)
    expect_equal(ans, 0)
})
