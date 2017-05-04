source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of the optim() function in NIMBLE code")

example <- list(
    par = 1.234,
    fn = function(par) { 2.345 + sum((par - 3.456)^2) },
    nimFnRun = function(par = double(1)) {
        return(2.345 + sum((par - 3.456)^2))
        returnType(double(0))
    }
)

normalizeWhitespace <- function(lines) {
    line <- paste(lines, collapse = ' ')
    line <- gsub('\\s+', ' ', line)  # Shrink internal whitespace.
    line <- gsub('^\\s+', '', line)  # Remove leading whitespace.
    line <- gsub('\\s+$', '', line)  # Remove trailing whitespace.
    return(line)
}

test_that("normalizeWhiteSpace() works", {
    expect_equal(normalizeWhitespace('  a b   cde  fg hi'), 'a b cde fg hi')
})

test_that("fakeOptim() behaves as expected", {
    actual <- fakeOptim(example$par, example$fn)
    expected <- optimResultNimbleList$new()
    expected$par <- example$par
    expected$value <- example$fn(expected$par)
    expect_equal(actual, expected)
})

test_that("nimbleFunction() replaces fakeOptim() with nimFakeOptim()", {
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(fakeOptim(par))
            returnType(optimResultNimbleList())
        }
    )()
    expect_equal(normalizeWhitespace(deparse(nimFun$run@.Data)), 'function (par) { return(nimFakeOptim(par)) }')
})

test_that("nimbleFunction() with optim() evaluates correctly", {
    fun <- function(par) {
        return(fakeOptim(par, example$fn))
    }
    nimFn <- nimbleFunction(run = example$nimFnRun)
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(fakeOptim(par, nimFn))
            returnType(optimResultNimbleList())
        }
    )()
    expect_equal(nimFun$run(example$par), fun(example$par))
})

test_that("compileNimble of optim stub works", {
    fun <- function(par) {
        return(fakeOptim(par, example$fn))
    }
    nimFn <- nimbleFunction(run = example$nimFnRun)
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(fakeOptim(par, nimFn))
            returnType(optimResultNimbleList())
        }
    )()
    compiledFun <- compileNimble(nimFun, showCompilerOutput = TRUE, dirName = '~/tmp')
    expect_equal(nimFun$run(example$par), compiledFun$run(example$par))
})
