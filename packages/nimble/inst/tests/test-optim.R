source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of the optim() function in NIMBLE code")

# Test helper to verify code.
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

example <- list(
    par = 1.234,
    fn = function(par) { 2.345 + sum((par - 3.456)^2) },
    nimFnRun = function(par = double(1)) {
        return(2.345 + sum((par - 3.456)^2))
        returnType(double(0))
    }
)

test_that("nimbleFunction() replaces optim() with nimOptim()", {
    nimFn <- nimbleFunction(run = example$nimFnRun)
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(optim(par, nimFn))
            returnType(optimResultNimbleList())
        }
    )()
    expect_equal(normalizeWhitespace(deparse(nimFun$run@.Data)),
                 'function (par) { return(nimOptim(par, nimFn)) }')
})

test_that("nimbleFunction with optim() agrees with R behavior", {
    fun <- function(par) {
        return(optim(par, example$fn))
    }
    nimFn <- nimbleFunction(run = example$nimFnRun)
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(optim(par, nimFn))
            returnType(optimResultNimbleList())
        }
    )()
    expect_equal(nimFun$run(example$par), fun(example$par))
})

test_that("compiled function with optim() agrees with nimbleFunction behavior", {
    nimFn <- nimbleFunction(run = example$nimFnRun)
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(optim(par, nimFn))
            returnType(optimResultNimbleList())
        }
    )()
    compiledFun <- compileNimble(nimFun, showCompilerOutput = TRUE)
    expect_failure(
        expect_equal(nimFun$run(example$par), compiledFun$run(example$par))
    )
})
