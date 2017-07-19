## This runs benchmarks of compiled NIMBLE code. It is located in tests/ so
## that it can continue to be tested, but it most useful to run this manually
## via `make benchmark` or `make profile`.
##
## Environment Variables:
##   NIMBLE_BENCHMARK_DIR - Directory where generated code is written.
##   NIMBLE_BENCHMARK_SEC - Minumum time to run each benchmark.

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context('Benchmarking Nimble code')
cat('\n')

## Computes number of iterations per second.
## This increase iterations until test takes around 1 sec.
benchmarkIters <- function(fun, minIters = 100L) {
    minRunTimeSec = 0.1
    if (nchar(Sys.getenv('NIMBLE_BENCHMARK_SEC'))) {
        minRunTimeSec <- as.double(Sys.getenv('NIMBLE_BENCHMARK_SEC'))
    }

    iters <- minIters
    elapsed <- fun(iters)
    while (elapsed < minRunTimeSec) {
        iters <- as.integer(iters * max(2, (minRunTimeSec / (elapsed + 1e-6)) ^ 0.8))
        elapsed <- fun(iters)
    }
    return(iters / elapsed)
}

## This makes DSL, C++, and possibly Tensorflow versions of a function.
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
    if (require('tensorflow')) {
        versions$tf <- withTensorflowEnabled(
            compileNimble(
                nimbleFunction(run = run, name = paste0('tf_', name)),
                projectName = 'tf', dirName = dirName))
    }
    return(versions)
}

matrixSizes <- c(1, 2, 3, 4, 5, 8, 10, 16, 32, 64, 100, 128, 256)
vectorSizes <- c(1, 2, 3, 4, 5, 8, 10, 16, 32, 64, 100, 128, 256,
                 512, 1000, 1024, 2048, 4096, 8192, 10000)

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

    cat('-------------------------------------------------\n')
    cat('Benchmarking Matrix Arithmetic\n')
    cat('  K   M   N DSL calls/ms C++ calls/ms TF calls/ms\n')
    for (K in matrixSizes) {
        M <- K
        N <- K
        x <- matrix(rnorm(K * M), K, M)
        y <- matrix(rnorm(M * N), M, N)
        nimPerSec <- benchmarkIters(function(iters) versions$dsl(x, y, iters))
        cPerSec <- benchmarkIters(function(iters) versions$cpp(x, y, iters))
        tfPerSec <- if (is.null(versions$tf)) NA else {
            benchmarkIters(function(iters) versions$tf(x, y, iters))
        }
        cat(sprintf('%3d %3d %3d %12.1f %12.1f %11.1f\n',
                    K, M, N, nimPerSec / 1e3, cPerSec / 1e3, tfPerSec / 1e3))
    }
    cat('-------------------------------------------------\n')
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
    cat('  K   M   N DSL calls/ms C++ calls/ms TF calls/ms\n')
    for (K in matrixSizes) {
        M <- K
        N <- K
        x <- matrix(rnorm(K * M), K, M)
        y <- matrix(rnorm(M * N), M, N)
        nimPerSec <- benchmarkIters(function(iters) versions$dsl(x, y, iters))
        cPerSec <- benchmarkIters(function(iters) versions$cpp(x, y, iters))
        tfPerSec <- if (is.null(versions$tf)) NA else {
            benchmarkIters(function(iters) versions$tf(x, y, iters))
        }
        cat(sprintf('%3d %3d %3d %12.1f %12.1f %11.1f\n',
                    K, M, N, nimPerSec / 1e3, cPerSec / 1e3, tfPerSec / 1e3))
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

    cat('-------------------------------------------\n')
    cat('Benchmarking Special Functions\n')
    cat('    N DSL calls/ms C++ calls/ms TF calls/ms\n')
    for (N in vectorSizes) {
        x <- exp(-rnorm(N))
        nimPerSec <- benchmarkIters(function(iters) versions$dsl(x, iters))
        cPerSec <- benchmarkIters(function(iters) versions$cpp(x, iters))
        tfPerSec <- if (is.null(versions$tf)) NA else {
            benchmarkIters(function(iters) versions$tf(x, iters))
        }
        cat(sprintf('%5d %12.1f %12.1f %11.1f\n',
                    N, nimPerSec / 1e3, cPerSec / 1e3, tfPerSec / 1e3))
    }
    cat('-------------------------------------------\n')
})
