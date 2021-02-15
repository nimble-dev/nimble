source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = FALSE)


context("Testing of the optim() function in NIMBLE code")

# The methods "SANN" and "Brent" are not supported.
methodsAllowingGradient <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B")
methodsAllowingBounds <- c("L-BFGS-B")

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

# This helper checks for agreement of the results to optim() and nimOptim() modulo some type mismatch.
expect_r_and_dsl_agree <- function(result_r, result_dsl, info = NULL) {
    # The return types differ.
    expect_equal(class(result_r), 'list', info = info)
    expect_equal(as.character(class(result_dsl)), 'OptimResultNimbleList', info = info)
    # But the values mostly agree.
    expect_equal(result_r$par, result_dsl$par, info = info)
    expect_equal(result_r$convergence, result_dsl$convergence, info = info)
    expect_equal(result_r$value, result_dsl$value, info = info)
    expect_equal(unname(result_r$counts), result_dsl$counts, info = info)
}

test_that("nimOptim() behaves mostly like optim()", {
    par <- c(1, 2, 3, 4)
    fn <- function(x) { sum(x ^ 2) }
    result_r <- optim(par, fn)
    result_dsl <- nimOptim(par, fn)
    expect_r_and_dsl_agree(result_r, result_dsl)
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
    expect_true(is(nimCaller$run(par), "OptimResultNimbleList"))
})

test_that("when an RCfunction optim()izes an RCfunction, the R and DSL behavior mostly agree", {
    # Define R versions.
    callee <- function(par) { return(sum(par ^ 2)) }
    caller <- function(par, method) { return(optim(par, callee, method = method)) }
    # Define DSL versions.
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimCallee, method = method))
            returnType(optimResultNimbleList())
        }
    )
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test agreement.
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_r <- caller(par, method)
        result_dsl <- nimCaller(par, method)
        expect_r_and_dsl_agree(result_r, result_dsl, info = paste(' where method =', method))
    }
    expect_error(caller(par, "bogus-method"))
    expect_error(nimCaller(par, "bogus-method"))
})

test_that("when a nimbleFunction optim()izes an RCfunction, the R and DSL behavior mostly agree", {
    # Define R versions.
    callee <- function(par) { return(sum(par ^ 2)) }
    caller <- function(par, method) { return(optim(par, callee, method = method)) }
    # Define DSL versions.
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        setup = TRUE,
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimCallee, method = method))
            returnType(optimResultNimbleList())
        }
    )()
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test agreement.
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_r <- caller(par, method)
        result_dsl <- nimCaller$run(par, method)
        expect_r_and_dsl_agree(result_r, result_dsl, info = paste(' where method =', method))
    }
    expect_error(caller(par, "bogus-method"))
    expect_error(nimCaller(par, "bogus-method"))
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
    # Test agreement.
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_r <- caller(par, method)
        result_dsl <- nimCaller$run(par, method)
        expect_r_and_dsl_agree(result_r, result_dsl, info = paste(' where method =', method))
    }
    expect_error(caller(par, "bogus-method"))
    expect_error(nimCaller$run(par, "bogus-method"))
})

test_that("when a nimbleFunction optim()izes a nimbleFunction, the R and DSL behavior mostly agree", {
    # Define R versions.
    callee <- function(par) { return(sum(par ^ 2)) }
    caller <- function(par, method) { return(optim(par, callee, method = method)) }
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
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimCallee$run, method = method))
            returnType(optimResultNimbleList())
        }
    )()
    # Test agreement.
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_r <- caller(par, method)
        result_dsl <- nimCaller$run(par, method)
        expect_r_and_dsl_agree(result_r, result_dsl, info = paste(' where method =', method))
    }
    expect_error(caller(par, "bogus-method"))
    expect_error(nimCaller$run(par, "bogus-method"))
})

test_that("when an RCfunction optim()izes an RCfunction, the DSL and C++ behavior agree", {
    nimCallee <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par ^ 2))
            returnType(double(0))
        }
    )
    nimCaller <- nimbleFunction(
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimCallee, method = method))
            returnType(optimResultNimbleList())
        }
    )
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test agreement.
    compiledCaller <- compileNimble(nimCaller)
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_dsl <- nimCaller(par, method)
        result_cpp <- compiledCaller(par, method)
        expect_equal(result_dsl, result_cpp, info = paste(' where method =', method))
    }
    expect_error(nimCaller(par, "bogus-method"))
    expect_error(compiledCaller(par, "bogus-method"))
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
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimCallee, method = method))
            returnType(optimResultNimbleList())
        }
    )()
    temporarilyAssignInGlobalEnv(nimCallee)  # Work around scoping issues.
    # Test agreement.
    compiledCaller <- compileNimble(nimCaller)
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_dsl <- nimCaller$run(par, method)
        result_cpp <- compiledCaller$run(par, method)
        expect_equal(result_dsl, result_cpp, info = paste(' where method =', method))
    }
    expect_error(nimCaller$run(par, "bogus-method"))
    expect_error(compiledCaller$run(par, "bogus-method"))
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
        run = function(par = double(1), method = character(0)) {
            return(optim(par, nimCallee$run, method = method))
            returnType(optimResultNimbleList())
        }
    )()
    # Test agreement.
    compiledCaller <- compileNimble(nimCaller)
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_dsl <- nimCaller$run(par, method)
        result_cpp <- compiledCaller$run(par, method)
        expect_equal(result_dsl, result_cpp, info = paste(' where method =', method))
    }
    expect_error(nimCaller$run(par, "bogus-method"))
    expect_error(compiledCaller$run(par, "bogus-method"))
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
    compiledCaller <- compileNimble(nimCaller)
    # Test agreement.
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_dsl <- nimCaller(par, method)
        result_cpp <- compiledCaller(par, method)
        expect_equal(result_dsl, result_cpp, info = paste(' where method =', method))
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
    compiledCaller <- compileNimble(nimCaller)
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_dsl <- nimCaller$run(par, method)
        result_cpp <- compiledCaller$run(par, method)
        expect_equal(result_dsl, result_cpp, info = paste(' where method =', method))
    }
    expect_error(nimCaller$run(par, "bogus-method"))
    expect_error(compiledCaller$run(par, "bogus-method"))
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
    compiledCaller <- compileNimble(nimCaller)
    # Test approximate agreement (i.e. that most fields agree).
    par <- c(1.2, 3.4)
    for (method in methodsAllowingGradient) {
        result_dsl <- nimCaller$run(par, method)
        result_cpp <- compiledCaller$run(par, method)
        expect_equal(result_dsl, result_cpp, info = paste(' where method =', method))
    }
    expect_error(nimCaller$run(par, "bogus-method"))
    expect_error(compiledCaller$run(par, "bogus-method"))
})

test_that("when optim() is called with bounds, behavior is the same among R, DSL, and C++", {
    # Only method L-BFGS-B supports bounds.
    # Define R functions.
    fn <- function(par) { sum((par - c(5, 7)) ^ 2) }
    optimizer <- function(lower, upper) {
        return(optim(c(1, 2), fn, method = 'L-BFGS-B', lower = lower, upper = upper))
    }
    # Define DSL functions.
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum((par - c(5, 7)) ^ 2))
            returnType(double(0))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)
    nimOptimizer <- nimbleFunction(
        run = function(lower = double(1), upper = double(1)) {
            par <- c(1, 2)  # FIXME compilation fails when this is inline.
            return(optim(par, nimFn, method = 'L-BFGS-B', lower = lower, upper = upper))
            returnType(optimResultNimbleList())
        }
    )
    compiledOptimizer <- compileNimble(nimOptimizer)
    
    expect_agreement <- function(lower, upper) {
        info <- paste(' where lower =', capture.output(dput(lower, control = c())),
                      ', upper =', capture.output(dput(upper, control = c())))
        result_r <- optimizer(lower, upper)
        result_dsl <- nimOptimizer(lower, upper)
        result_cpp <- compiledOptimizer(lower, upper)
        expect_r_and_dsl_agree(result_r, result_dsl, info = info)
        expect_equal(result_dsl, result_cpp, info = info)
    }
    
    # Test with many configurations.
    expect_agreement(lower = -Inf, upper = Inf)
    expect_agreement(lower = -Inf, upper = 2)
    expect_agreement(lower = 6, upper = Inf)
    expect_agreement(lower = c(6, 6), upper = Inf)
    if(RUN_FAILING_TESTS) {  # TODO Fix these failing tests.
        expect_agreement(lower = 2, upper = 3)  # FAILS dsl ! = cpp only in .count field
        expect_agreement(lower = c(1, 2), upper = c(3, 4))  # FAILS dsl != cpp only in .count field
        expect_agreement(lower = -Inf, upper = c(6, 6))  # FAILS dsl != cpp only in .count field
    }
})

test_that("when optim() is called with bounds and gradient, behavior is the same among R, DSL, and C++", {
    # Only method L-BFGS-B supports bounds.
    # Define R functions.
    fn <- function(par) { sum((par - c(5, 7)) ^ 2) }
    gr <- function(par) { 2 * (par - c(5, 7)) }
    optimizer <- function(lower, upper) {
        return(optim(c(1, 2), fn, gr, method = 'L-BFGS-B', lower = lower, upper = upper))
    }
    # Define DSL functions.
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum((par - c(5, 7)) ^ 2))
            returnType(double(0))
        }
    )
    nimGr <- nimbleFunction(
        run = function(par = double(1)) {
            return(2 * (par - c(5, 7)))
            returnType(double(1))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)
    temporarilyAssignInGlobalEnv(nimGr)
    nimOptimizer <- nimbleFunction(
        run = function(lower = double(1), upper = double(1)) {
            par <- c(1, 2)  # FIXME compilation fails when this is inline.
            return(optim(par, nimFn, nimGr, method = 'L-BFGS-B', lower = lower, upper = upper))
            returnType(optimResultNimbleList())
        }
    )
    compiledOptimizer <- compileNimble(nimOptimizer)
    
    expect_agreement <- function(lower, upper) {
        info <- paste(' where lower =', capture.output(dput(lower, control = c())),
                      ', upper =', capture.output(dput(upper, control = c())))
        result_r <- optimizer(lower, upper)
        result_dsl <- nimOptimizer(lower, upper)
        result_cpp <- compiledOptimizer(lower, upper)
        expect_r_and_dsl_agree(result_r, result_dsl, info = info)
        expect_equal(result_dsl, result_cpp, info = info)
    }
    
    # Test with many configurations.
    expect_agreement(lower = -Inf, upper = Inf)
    expect_agreement(lower = -Inf, upper = 2)
    expect_agreement(lower = 6, upper = Inf)
    expect_agreement(lower = c(6, 6), upper = Inf)
    if(RUN_FAILING_TESTS) {  # TODO Fix these failing tests.
        expect_agreement(lower = 2, upper = 3)  # FAILS dsl ! = cpp only in .count field
        expect_agreement(lower = c(1, 2), upper = c(3, 4))  # FAILS dsl != cpp only in .count field
        expect_agreement(lower = -Inf, upper = c(6, 6))  # FAILS dsl != cpp only in .count field
    }
})

test_that("when optim() is called with control, behavior is the same among R, DSL, and C++", {
    # Define R functions.
    fn <- function(par) { sum((par - c(5, 7)) ^ 2) }
    gr <- function(par) { 2 * (par - c(5, 7)) }
    optimizer <- function(method, maxit, type) {
        control <- list(maxit = maxit, type = type)
        return(optim(c(1, 2), fn, gr, method = method, control = control))
    }
    # Define DSL functions.
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum((par - c(5, 7)) ^ 2))
            returnType(double(0))
        }
    )
    nimGr <- nimbleFunction(
        run = function(par = double(1)) {
            return(2 * (par - c(5, 7)))
            returnType(double(1))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)
    temporarilyAssignInGlobalEnv(nimGr)
    nimOptimizer <- nimbleFunction(
        run = function(method = character(0), maxit = integer(0), type = integer(0)) {
            control <- optimDefaultControl()
            control$maxit <- maxit
            control$type <- type
            par <- c(1, 2)  # FIXME compilation fails when this is inline.
            return(optim(par, nimFn, nimGr, method = method, control = control))
            returnType(optimResultNimbleList())
        }
    )
    compiledOptimizer <- compileNimble(nimOptimizer)

    expect_agreement <- function(method, ...) {
        control_nondefault <- list(...)
        control <- list(maxit = 100, type = 1)  # Defaults.
        for (name in names(control_nondefault)) {
            control[[name]] = control_nondefault[[name]]
        }
        result_r <- optimizer(method, control$maxit, control$type)
        result_dsl <- nimOptimizer(method, control$maxit, control$type)
        result_cpp <- compiledOptimizer(method, control$maxit, control$type)
        info <- capture.output(dput(list(method = method, control = control), control = c()))
        expect_r_and_dsl_agree(result_r, result_dsl, info = info)
        expect_equal(result_dsl, result_cpp, info = info)
    }

    # Test with many configurations.
    expect_agreement(method = "Nelder-Mead", maxit = 5)
    expect_agreement(method = "BFGS", maxit = 5)
    expect_agreement(method = "CG", maxit = 5)
    expect_agreement(method = "L-BFGS-B", maxit = 5)
    expect_agreement(method = "CG", type = 1)  # Fletcher-Reeves
    expect_agreement(method = "CG", type = 2)  # Polak-Ribiere
    expect_agreement(method = "CG", type = 3)  # Beal-Sorenson
})

test_that("optim() minimizes with fnscale = 1 and maximizes when fnscale = -1", {
    # Define R functions with exactly one global minimumm and one global maximum.
    fn <- function(par) {
        return(sum(par) * exp(-sum(par ^ 2) / 2))
    }
    optimizer <- function(method, fnscale) {
        control <- list(fnscale = fnscale)
        return(optim(c(0.1, -0.1), fn, method = method, control = control))
    }
    # Define DSL functions.
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par) * exp(-sum(par ^ 2) / 2))
            returnType(double(0))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)
    nimOptimizer <- nimbleFunction(
        run = function(method = character(0), fnscale = double(0)) {
            control <- optimDefaultControl()
            control$fnscale <- fnscale
            par <- c(0.1, -0.1)  # FIXME compilation fails when this is inline.
            return(optim(par, nimFn, method = method, control = control))
            returnType(optimResultNimbleList())
        }
    )
    compiledOptimizer <- compileNimble(nimOptimizer)
    # Test with many methods and fnscales.
    for (method in methodsAllowingGradient) {
        for (fnscale in c(-1, 1)) {
            result_r <- optimizer(method, fnscale)
            result_dsl <- nimOptimizer(method, fnscale)
            result_cpp <- compiledOptimizer(method, fnscale)
            info <- paste(' where method =', method, ', fnscale =', fnscale)
            expect_r_and_dsl_agree(result_r, result_dsl, info = info)
            expect_equal(result_dsl, result_cpp, info = info)
        }
    }
})

test_that("optim() with gradient minimizes with fnscale = 1 and maximizes when fnscale = -1", {
    # Define R functions with exactly one global minimumm and one global maximum.
    fn <- function(par) {
        return(sum(par) * exp(-sum(par ^ 2) / 2))
    }
    gr <- function(par) {
        return(sum(par) * exp(-sum(par ^ 2) / 2) * (rep(1, length(par)) - par))
    }
    optimizer <- function(method, fnscale) {
        control <- list(fnscale = fnscale)
        return(optim(c(0.1, -0.1), fn, gr, method = method, control = control))
    }
    # Define DSL functions.
    nimFn <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par) * exp(-sum(par ^ 2) / 2))
            returnType(double(0))
        }
    )
    nimGr <- nimbleFunction(
        run = function(par = double(1)) {
            return(sum(par) * exp(-sum(par ^ 2) / 2) * (rep(1, length(par)) - par))
            returnType(double(1))
        }
    )
    temporarilyAssignInGlobalEnv(nimFn)
    temporarilyAssignInGlobalEnv(nimGr)
    nimOptimizer <- nimbleFunction(
        run = function(method = character(0), fnscale = double(0)) {
            control <- optimDefaultControl()
            control$fnscale <- fnscale
            par <- c(0.1, -0.1)  # FIXME compilation fails when this is inline.
            return(optim(par, nimFn, nimGr, method = method, control = control))
            returnType(optimResultNimbleList())
        }
    )
    compiledOptimizer <- compileNimble(nimOptimizer)
    # Test with many methods and fnscales.
    for (method in methodsAllowingGradient) {
        for (fnscale in c(-1, 1)) {
            result_r <- optimizer(method, fnscale)
            result_dsl <- nimOptimizer(method, fnscale)
            result_cpp <- compiledOptimizer(method, fnscale)
            info <- paste(' where method =', method, ', fnscale =', fnscale)
            expect_r_and_dsl_agree(result_r, result_dsl, info = info)
            expect_equal(result_dsl, result_cpp, info = info)
        }
    }
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)