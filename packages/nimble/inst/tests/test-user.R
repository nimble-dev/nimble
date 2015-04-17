#source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))


# tmp
library(testthat)
context("Testing of user-supplied distributions and functions in BUGS code")

## User-supplied distributions

# tmp
library(nimble, lib.loc = '/tmp/nim41')

dbl <- nimbleFunction(
    run = function(x = double(0)) {
        returnType(double(0))
        return(2*x)
    }
    )

## vector input and output
vecdbl <- nimbleFunction(
    run = function(x = double(1)) {
        returnType(double(1))
        return(2*x)
    }
    )

## two arguments to the nimbleFunction
dblSum <- nimbleFunction(
    run = function(x = double(0), y = double(0)) {
        returnType(double(0))
        return(2*(x+y))
    }
    )


code <- nimbleCode({
    x ~ dnorm(0, 1)
    dx ~ dnorm(dbl(x), .01)
    y[1:K] ~ dmnorm(vecdbl(mu[1:K]), cov = .01*I[1:K, 1:K])
    mu[1:K] ~ dmnorm(zeros[1:K], cov = I[1:K, 1:K])
    z ~ dnorm(0, 1)
    dz ~ dnorm(dblSum(x, z), .01)
    # vectorized fun applied to scalar nodes-based variable

    for(i in 1:K) {
        theta[i] ~ dnorm(0, 1)
    }
    w[1:K] ~ dmnorm(vecdbl(theta[1:K]), cov = .01*I[1:K, 1:K])
})

K <- 3
m <- nimbleModel(code, inits = list(x = 0.25, y = 1:K, mu = 1:K,
                           z = 0.5, theta = rep(.5, K), w = rep(1, K)),
                 constants = list(K = K, zeros = rep(0, K), I = diag(K)))

cm <- compileNimble(m)

if(F) {
set.seed(0)
simulate(m)
set.seed(0)
simulate(cm)

for(var in c('dx', 'y', 'dz', 'w')) {
    try(test_that("Test that R and C models agree with user-suppled functions: ",
                  expect <- that(get(var, m), equals(get(var, cm),
                                             info = paste0(var, " values differ")))))
}

}
                   
    

    
