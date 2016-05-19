## File for testing distributions provided by NIMBLE

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
context('Testing NIMBLE distributions')

## mvt

x <- c(1, 1, 2)
mn <- c(1, 2, 3)
cov <- diag(c(1, 2, 3))
cov[1, 3] <- cov[3, 1] <- 0.1
df <- 5


## test use directly from R

truth <- mvtnorm::dmvt(x, delta = mn, sigma = cov, df = df, log = FALSE)

try(test_that("dmvt_chol calculates density correctly in R: ",
              expect_equal(dmvt_chol(x, mn, chol(cov), df, prec_param = FALSE),
                          truth,
                          info = paste0("incorrect dmvt calculation in R"))))

## test use through nimble function

nf <- nimbleFunction(
    run = function(x = double(1), mn = double(1),
                   cov = double(2), df = double(0)) {
        returnType(double(0))
        ch <- chol(cov)
        out <- dmvt_chol(x = x, mean = mn, cholesky = ch,
                         df = df, prec_param = FALSE, log = FALSE)
        return(out)
    }
)

cnf <- compileNimble(nf)

try(test_that("Test that dmvt_chol works correctly in R nimble function: ",
                  expect_equal(nf(x, mn, cov, df), (truth), 
                  info = paste0("incorrect dmvt value in R nimble function"))))

try(test_that("Test that dmvt_chol works correctly in compiled nimble function: ",
                  expect_equal(cnf(x, mn, cov, df), (truth), 
                  info = paste0("incorrect dmvt value in compiled nimble function"))))

## test use in model

mvt_code <- nimbleCode({
    x[1:3] ~ dmvt(mn[], cov = cov[,], df = df)
})

mvt_model <- nimbleModel(mvt_code, constants = list(mn = mn, cov = cov, prec = FALSE, df = df))

mvt_model$x <- x

try(test_that("Test that dmvt calculation is correct in model likelihood calculation: ",
              expect_equal(exp(mvt_model$calculate()), (truth),
                          info = paste0("incorrect likelihood value for dmvt"))))

c_mvt_model <- compileNimble(mvt_model)
c_mvt_model$x
try(test_that("Test that dmvt (compiled) calculation is correct in model likelihood calculation: ",
              expect_equal(exp(c_mvt_model$calculate()), (truth),
                          info = paste0("incorrect likelihood value for dmvt (compiled)"))))

## random sampling
reference_samps <- mvtnorm::rmvt(n = 10000, delta = mn, sigma = cov, df = df)
r_samps <- t(replicate(10000, rmvt_chol(n = 1, mn, chol(cov), df, prec_param = FALSE)))

try(test_that("Test that random samples (R) have correct mean: ",
              expect_equal(colMeans(r_samps), (colMeans(reference_samps)),
                                                    tol = 0.01,
                          info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (R) have correct covariance: ",
              expect_equal(cov(r_samps), (cov(reference_samps)),
                                               tol = 0.1,
                          info = "Difference in covs exceeds tolerance")))


nf_sampling <- nimbleFunction(
    run = function(mn = double(1), cov = double(2), df = double(0)) {
        returnType(double(1))
        ch <- chol(cov)
        out <- rmvt_chol(n = 1, mean = mn, cholesky = ch,
                         df = df, prec_param = FALSE)
        return(out)
    }
)

nf_samps <- t(replicate(10000, nf_sampling(mn, cov, df)))

try(test_that("Test that random samples (nf) have correct mean: ",
              expect_equal(colMeans(nf_samps), (colMeans(reference_samps)),
                                                    tol = 0.01,
                          info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (nf) have correct covariance: ",
              expect_equal(cov(nf_samps), (cov(reference_samps)),
                                                tol = 0.1,
                          info = "Difference in covs exceeds tolerance")))

## sampling via `simulate`
simul_samp <- function(model) {
    model$simulate()
    return(model$x)
}

## this is rather slow, so we'll have to test with fewer samples
simul_samps <- t(replicate(1000, simul_samp(c_mvt_model)))

try(test_that("Test that random samples (simulate) have correct mean: ",
              expect_equal(colMeans(simul_samps), (colMeans(reference_samps)),
                                                    tol = 0.1,
                          info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (simulate) have correct covariance: ",
              expect_equal(cov(simul_samps), (cov(reference_samps)),
                                                tol = 0.5,
                          info = "Difference in covs exceeds tolerance")))
