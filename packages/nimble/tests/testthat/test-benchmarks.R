## This runs benchmarks of compiled NIMBLE code. It is located in tests/ so
## that it can continue to be tested, but it most useful to run this manually
## via `make benchmark` or `make profile`.
##
## Environment Variables:
##   NIMBLE_BENCHMARK_DIR - Directory where generated code is written.
##   NIMBLE_BENCHMARK_SEC - Minumum time to run each benchmark.

source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

# This makes it easier to debug compiler and linker errors on travis.
if (nchar(Sys.getenv('CI'))) nimbleOptions(showCompilerOutput = TRUE)

RwarnLevel <- options('warn')$warn
options(warn = 1)

context('Benchmarking NIMBLE code')
cat('\n')

## Computes number of iterations per second.
## This increases iterations until test takes around 1 sec.
benchmarkIters <- function(fun, minIters = 10L) {
    minRunTimeSec = 0.1
    if (nchar(Sys.getenv('NIMBLE_BENCHMARK_SEC'))) {
        minRunTimeSec <- as.double(Sys.getenv('NIMBLE_BENCHMARK_SEC'))
    }

    # Check availability and warm-up.
    if (is.na(fun(1))) return(NA)

    iters <- minIters
    elapsed <- fun(iters)
    while (elapsed < minRunTimeSec && iters < 1000000) {
        iters <- as.integer(iters * max(3, min(100, (minRunTimeSec / (elapsed + 1e-6)) ^ 0.8)))
        elapsed <- fun(iters)
    }
    return(iters / elapsed)
}

## This makes DSL and C++ versions of a function.
## We give each a distinct name so that profiling tools like callgrind can disambiguate function symbols.
makeVersions <- function(name, run) {
    dirName <- NULL
    if (nchar(Sys.getenv('NIMBLE_BENCHMARK_DIR'))) {
        dirName <- Sys.getenv('NIMBLE_BENCHMARK_DIR')
    }

    versions <- list()
    versions$dsl <- nimbleFunction(run = run, name = paste0('dsl_', name))
    versions$cpp <- compileNimble(
        nimbleFunction(run = run, name = paste0('cpp_', name)),
        projectName = 'cpp', dirName = dirName)
    return(versions)
}

matrixSizes <- c(1, 10, 100, 1000)
vectorSizes <- c(1, 10, 100, 1000, 10000, 100000)

test_that('Benchmarking matrix arithmetic', {
    versions <- makeVersions(
        name = 'arithmetic',
        run = function(x = double(2), y = double(2), numIters = integer(0)) {
            result <- run.time({
                for (i in 1:numIters) {
                    z <- 2 * x + y
                    w <- x ^ 2 - y * z
                }
            })
            return(result)
            returnType(double(0))
        })

    cat('--------------------------------------------\n')
    cat('Benchmarking Matrix Arithmetic\n')
    cat('   M    N DSL ops/sec C++ ops/sec\n')
    for (M in matrixSizes) {
        N <- M
        x <- matrix(rnorm(M * N), M, N)
        y <- matrix(rnorm(M * N), M, N)
        nimPerSec <- benchmarkIters(function(iters) versions$dsl(x, y, iters))
        cPerSec <- benchmarkIters(function(iters) versions$cpp(x, y, iters))
        cat(sprintf('%4d %4d %11.2g %11.2g\n',
                    M, N, nimPerSec, cPerSec))
    }
    cat('--------------------------------------------\n')
})

test_that('Benchmarking matrix multiplication', {
    versions <- makeVersions(
        name = 'matmul',
        run = function(x = double(2), y = double(2), numIters = integer(0)) {
            result <- run.time({
                for (i in 1:numIters) {
                    z <- x %*% y
                }
            })
            return(result)
            returnType(double(0))
        })

    cat('-------------------------------------------------\n')
    cat('Benchmarking Matrix Multiplication\n')
    cat('   K    M    N DSL ops/sec C++ ops/sec/sec\n')
    for (K in matrixSizes) {
        M <- K
        N <- K
        x <- matrix(rnorm(K * M), K, M)
        y <- matrix(rnorm(M * N), M, N)
        nimPerSec <- benchmarkIters(function(iters) versions$dsl(x, y, iters))
        cPerSec <- benchmarkIters(function(iters) versions$cpp(x, y, iters))
        cat(sprintf('%4d %4d %4d %11.2g %11.2g\n',
                    K, M, N, nimPerSec, cPerSec))
    }
    cat('-------------------------------------------------\n')
})

test_that('Benchmarking vectorized special functions', {
    versions <- makeVersions(
        name = 'special',
        run = function(x = double(1), numIters = integer(0)) {
            result <- run.time({
                for (i in 1:numIters) {
                    y <- sqrt(x) + log(x) + exp(x) + cos(x) + sin(x) + lgamma(x) + expit(x)
                }
            })
            return(result)
            returnType(double(0))
        })

    cat('-----------------------------------------\n')
    cat('Benchmarking Special Functions\n')
    cat('     N DSL ops/sec C++ ops/sec\n')
    for (N in vectorSizes) {
        x <- exp(-rnorm(N))
        nimPerSec <- benchmarkIters(function(iters) versions$dsl(x, iters))
        cPerSec <- benchmarkIters(function(iters) versions$cpp(x, iters))
        cat(sprintf('%6d %11.1g %11.1g\n',
                    N, nimPerSec, cPerSec))
    }
    cat('-----------------------------------------\n')
})

options(warn = RwarnLevel)
