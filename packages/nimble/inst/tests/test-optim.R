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

test_that("nimbleFunction() replaces fakeOptim() with nimFakeOptim()", {
    nimFun <- nimbleFunction(
        run = function(x = double(0)) {
            return(fakeOptim(x))
            returnType(double(0))
        }
    )
    expect_equal(normalizeWhitespace(deparse(nimFun)), 'function (x) { return(nimFakeOptim(x)) }')
})

test_that("nimbleFunction() with optim() evaluates correctly", {
    fun <- function(x) {
        return(fakeOptim(x))
    }
    nimFun <- nimbleFunction(
        run = function(x = double(0)) {
            return(fakeOptim(x))
            returnType(double(0))
        }
    )
    ans <- nimFun(0)
    expect_equal(nimFun(0), 0)
})

test_that("compileNimble of optim stub works", {
    nimFun <- nimbleFunction(
        run = function(x = double(0)) {
            return(fakeOptim(x))
            returnType(double(0))
        }
    )
    compiledFun <- compileNimble(nimFun)
    ans <- compiledFun(0)
    expect_equal(ans, 0)
})
