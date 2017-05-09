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

test_that("nimOptim() behaves mostly like optim()", {
    par <- c(1, 2, 3, 4)
    fn <- function(x) { sum(x ^ 2) }
    expected <- optim(par, fn)
    actual <- nimOptim(par, fn)
    # The return types differ.
    expect_equal(class(expected), 'list')
    expect_equal(as.character(class(actual)), 'OptimResultNimbleList')
    # But the values mostly agree.
    expect_equal(actual$par, expected$par)
    expect_equal(actual$convergence, expected$convergence)
    expect_equal(actual$value, expected$value)
    expect_equal(length(actual$counts), length(expected$counts))
    expect_equal(actual$counts[1], expected$counts[[1]])  # Note the indexing disagreement.
})

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
    # Work around scoping issues.
    assign('nimCallee', nimCallee, envir = .GlobalEnv)
    on.exit(remove('nimCallee', envir = .GlobalEnv), add = TRUE)
    par <- c(1.2, 3.4)
    nimCaller$run(par)
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
    # Work around scoping issues.
    assign('nimCallee', nimCallee, envir = .GlobalEnv)
    on.exit(remove('nimCallee', envir = .GlobalEnv), add = TRUE)
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    expected <- caller(par)
    actual <- nimCaller$run(par)
    expect_equal(actual$par, expected$par)
    expect_equal(actual$convergence, expected$convergence)
    expect_equal(actual$value, expected$value)
    expect_equal(length(actual$counts), length(expected$counts))
    expect_equal(actual$counts[1], expected$counts[[1]])  # Note the indexing disagreement.
})

test_that("when a nimbleFunction optim()izes an RCfunction, the DSL and C++ behavior agree", {
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
    # Work around scoping issues.
    assign('nimCallee', nimCallee, envir = .GlobalEnv)
    on.exit(remove('nimCallee', envir = .GlobalEnv), add = TRUE)
    # Test agreement.
    par <- c(1.2, 3.4)
    compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
    expected <- nimCaller$run(par)
    actual <- compiledCaller$run(par)
    actual$counts[2] <- NA  # TODO remove this kludge by correctly treating NA in C++ -> R conversion.
    expect_equal(actual, expected)
})

test_that("when a nimbleFunction optim()izes a nimbleFunction, the DSL and C++ behavior agree", {
    nimCalleeGen <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    # Work around scoping issues.
    assign('nimCalleeGen', nimCalleeGen, envir = .GlobalEnv)
    on.exit(remove('nimCalleeGen', envir = .GlobalEnv), add = TRUE)
    nimCaller <- nimbleFunction(
        setup = function() {
            nimCallee <- nimCalleeGen()
        },
        run = function(par = double(1)) {
            return(optim(par, nimCallee$run))
            returnType(optimResultNimbleList())
        }
    )()
    # Test agreement.
    compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
    par <- c(1.2, 3.4)
    expected <- nimCaller$run(par)
    actual <- compiledCaller$run(par)
    actual$counts[2] <- NA  # TODO remove this kludge by correctly treating NA in C++ -> R conversion.
    expect_equal(actual, expected)
})
