source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of the optim() function in NIMBLE code")

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
    par <- 1.234
    actual <- fakeOptim(par)
    expected <- optimResultNimbleList$new()
    expected$par <- par
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
        return(fakeOptim(par))
    }
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(fakeOptim(par))
            returnType(optimResultNimbleList())
        }
    )()
    par <- 1.234
    expect_equal(nimFun$run(par), fun(par))
})

test_that("compileNimble of optim stub works", {
    nimFun <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(fakeOptim(par))
            returnType(optimResultNimbleList())
        }
    )()
    # compiledFun <- compileNimble(nimFun)
    compiledFun <- compileNimble(nimFun, showCompilerOutput = TRUE, dirName = '~/tmp')
    par <- 1.234
    expect_equal(nimFun$run(par), compiledFun$run(par))
})
