source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of the optim() function in NIMBLE code")

# The method "SANN" uses the gr arg to specify a transition kernel, so we test it independently.
methodsAllowingGradient <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")

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
    expect_equal(actual$counts, unname(expected$counts))
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
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    par <- c(1.2, 3.4)
    nimCaller$run(par)
})

test_that("when an RCfunction optim()izes an RCfunction, the R and DSL behavior mostly agree", {
    # Define R versions.
    callee <- function(par) { return(sum(par ^ 2)) }
    caller <- function(par) { return(optim(par, callee)) }
    # Define DSL versions.
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        run = function(par = double(1)) {
            return(optim(par, nimCallee))
            returnType(optimResultNimbleList())
        }
    )
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    expected <- caller(par)
    actual <- nimCaller(par)
    expect_equal(actual$par, expected$par)
    expect_equal(actual$convergence, expected$convergence)
    expect_equal(actual$value, expected$value)
    expect_equal(actual$counts, unname(expected$counts))
})

test_that("when a nimbleFunction optim()izes an RCfunction, the R and DSL behavior mostly agree", {
    # Define R versions.
    callee <- function(par) { return(sum(par ^ 2)) }
    caller <- function(par) { return(optim(par, callee)) }
    # Define DSL versions.
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
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    expected <- caller(par)
    actual <- nimCaller$run(par)
    expect_equal(actual$par, expected$par)
    expect_equal(actual$convergence, expected$convergence)
    expect_equal(actual$value, expected$value)
    expect_equal(actual$counts, unname(expected$counts))
})

test_that("when a nimbleFunction optim()izes an RCfunction with gradient, the R and DSL behavior mostly agree", {
    # Define R versions.
    fn <- function(par) { return(sum(par ^ 2)) }
    gr <- function(par) { return(2 * par) }
    caller <- function(par, method) {
        return(optim(par, fn, gr, method = method))
    }
    # Define DSL versions.
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimGr <- nimbleFunction(
        run = function(par = double(1)) {
            return(2 * par)
            returnType(double(1))
        }
    )
    nimCaller <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimFn, nimGr, method = method))
            returnType(optimResultNimbleList())
        }
    )()
    temporarilyAssignInGlobalEnv(nimFn)  # Work around scoping issues.
    temporarilyAssignInGlobalEnv(nimGr)  # Work around scoping issues.
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        info = paste(' where method =', method)
        expected <- caller(par, method)
        actual <- nimCaller$run(par, method)
        expect_equal(actual$par, expected$par, info = info)
        expect_equal(actual$convergence, expected$convergence, info = info)
        expect_equal(actual$value, expected$value, info = info)
        expect_equal(actual$counts, unname(expected$counts), info = info)
    }
    expect_error(caller(par, "bogus-method"))
    expect_error(nimCaller$run(par, "bogus-method"))
})

test_that("when a nimbleFunction optim()izes a nimbleFunction, the R and DSL behavior mostly agree", {
    # Define R versions.
    callee <- function(par) { return(sum(par ^ 2)) }
    caller <- function(par) { return(optim(par, callee)) }
    # Define DSL versions.
    nimCalleeGen <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    temporarilyAssignInGlobalEnv(nimCalleeGen)  # Work around scoping issues.
    nimCaller <- nimbleFunction(
        setup = function() {
            nimCallee <- nimCalleeGen()
        },
        run = function(par = double(1)) {
            return(optim(par, nimCallee$run))
            returnType(optimResultNimbleList())
        }
    )()
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    expected <- caller(par)
    actual <- nimCaller$run(par)
    expect_equal(actual$par, expected$par)
    expect_equal(actual$convergence, expected$convergence)
    expect_equal(actual$value, expected$value)
    expect_equal(actual$counts, unname(expected$counts))
})

test_that("when an RCfunction optim()izes an RCfunction, the DSL and C++ behavior agree", {
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        run = function(par = double(1)) {
            return(optim(par, nimCallee))
            returnType(optimResultNimbleList())
        }
    )
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test agreement.
    compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
    par <- c(1.2, 3.4)
    expected <- nimCaller(par)
    actual <- compiledCaller(par)
    expect_equal(actual, expected)
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
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test agreement.
    compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
    par <- c(1.2, 3.4)
    expected <- nimCaller$run(par)
    actual <- compiledCaller$run(par)
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
    temporarilyAssignInGlobalEnv(nimCalleeGen)  # Work around scoping issues.
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
    expect_equal(actual, expected)
})

test_that("when an RCfunction optim()izes an RCfunction with gradient, the DSL and C++ behavior mostly agree", {
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimGr <- nimbleFunction(
        run = function(par = double(1)) {
            return(2 * par)
            returnType(double(1))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)  # Work around scoping issues.
    temporarilyAssignInGlobalEnv(nimGr)  # Work around scoping issues.
    nimCaller <- nimbleFunction(
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimFn, nimGr, method = method))
            returnType(optimResultNimbleList())
        }
    )
    compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        info = paste(' where method =', method)
        expected <- nimCaller(par, method)
        actual <- compiledCaller(par, method)
        expect_equal(actual$par, expected$par, info = info)
        expect_equal(actual$convergence, expected$convergence, info = info)
        expect_equal(actual$value, expected$value, info = info)
        expect_equal(actual$counts, unname(expected$counts), info = info)
    }
    expect_error(nimCaller(par, "bogus-method"))
    expect_error(compiledCaller(par, "bogus-method"))
})

test_that("when a nimbleFunction optim()izes an RCfunction with gradient, the DSL and C++ behavior mostly agree", {
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimGr <- nimbleFunction(
        run = function(par = double(1)) {
            return(2 * par)
            returnType(double(1))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)  # Work around scoping issues.
    temporarilyAssignInGlobalEnv(nimGr)  # Work around scoping issues.
    nimCaller <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimFn, nimGr, method = method))
            returnType(optimResultNimbleList())
        }
    )()
    compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        info = paste(' where method =', method)
        expected <- nimCaller$run(par, method)
        actual <- compiledCaller$run(par, method)
        expect_equal(actual$par, expected$par, info = info)
        expect_equal(actual$convergence, expected$convergence, info = info)
        expect_equal(actual$value, expected$value, info = info)
        expect_equal(actual$counts, unname(expected$counts), info = info)
    }
    expect_error(nimCaller$run(par, "bogus-method"))
    expect_error(compiledCaller$run(par, "bogus-method"))
})

# FAILS, sometimes crashing rstudio
test_that("when an RCfunction optimizes an RCfunction with SANN, the DSL and C++ behavior mostly agree", {
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimGr <- nimbleFunction(
        run = function(par = double(1)) {
            # Work around https://github.com/nimble-dev/nimble/issues/351
            direction <- rnorm(length(par))
            return(direction)
            returnType(double(1))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)  # Work around scoping issues.
    temporarilyAssignInGlobalEnv(nimGr)  # Work around scoping issues.
    nimCaller <- nimbleFunction(
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimFn, nimGr, method = method))
            returnType(optimResultNimbleList())
        }
    )
    compiledCaller <- compileNimble(nimCaller, showCompilerOutput = TRUE)
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    expected <- nimCaller(par, method = "SANN")
    actual <- compiledCaller(par, method = "SANN")
    expect_equal(actual$par, expected$par)
    expect_equal(actual$convergence, expected$convergence)
    expect_equal(actual$value, expected$value)
    expect_equal(actual$counts, unname(expected$counts))
})
