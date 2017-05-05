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

# ----------------------------------------------------------------------------
# Simple tests that avoid R magic.

test_that("nimOptim() behaves mostly like optim()", {
    par <- c(1, 2, 3, 4)
    fn <- function(x) { sum(x ^ 2) }
    result <- optim(par, fn)
    nimResult <- nimOptim(par, fn)
    # The return types differ.
    expect_equal(class(result), 'list')
    expect_equal(as.character(class(nimResult)), 'OptimResultNimbleList')
    # But the values agree.
    for (n in names(result)) {
        expect_equal(result[[n]], nimResult[[n]])
    }
})

# We work around scoping issues with nimbleFunctions calling other nimbleFunctions by assigning nimCallee
# in the global scope. This is accomplished by the following assign() and below with_mock() calls.
assign('nimCallee', function() { fail('nimCalle should not have been called') }, envir = topenv())

test_that("nimbleFunction() replaces optim() with nimOptim()", {
    # Note that nimCallee may be undefined, since nimCaller$run() is never called.
    nimCaller <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(optim(par, nimCallee))
            returnType(optimResultNimbleList())
        }
    )()
    expect_equal(normalizeWhitespace(deparse(nimCaller$run@.Data)),
                 'function (par) { return(nimOptim(par, nimCallee)) }')
})

test_that("nimbleFunction with optim() runs", {
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(optim(par, nimCallee))
            returnType(optimResultNimbleList())
        }
    )()
    with_mock(nimCallee = nimCallee, {
        par <- c(1.2, 3.4)
        nimCaller$run(par)
    })
})

test_that("nimbleFunction with optim() mostly agrees with R behavior", {
    # Define R versions.
    callee <- function(par) { return(sum(par ^ 2)) }
    caller <- function(par) { return(optim(par, callee)) }
    # Define nimble versions.
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(optim(par, nimCallee))
            returnType(optimResultNimbleList())
        }
    )()
    # Test agreement.
    with_mock(nimCallee = nimCallee, {
        par <- c(1.2, 3.4)
        result <- caller(par)
        nimResult <- do.call(optimResultNimbleList$new, result)
        expect_equal(nimCaller$run(par), nimResult)
    })
})

test_that("nimbleFunction with optim() agrees with C++ behavior", {
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(optim(par, nimCallee))
            returnType(optimResultNimbleList())
        }
    )()
    # Test agreement.
    with_mock(nimCallee = nimCallee, {
        par <- c(1.2, 3.4)
        compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
        expect_equal(nimCaller$run(par), compiledCaller$run(par))
    })
})
