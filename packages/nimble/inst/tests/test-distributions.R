## File for testing distributions provided by NIMBLE

source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))
context('Testing NIMBLE distributions')

## mvt
x <- c(1, 1, 2)
mn <- c(1, 2, 3)
sc <- diag(c(1, 2, 3))
sc[1, 3] <- sc[3, 1] <- 0.1
df <- 5

## test use directly from R

truth <- mvtnorm::dmvt(x, delta = mn, sigma = sc, df = df, log = FALSE)

try(test_that("dmvt_chol calculates density correctly in R: ",
              expect_equal(dmvt_chol(x, mn, chol(sc), df, prec_param = FALSE),
                          truth,
                          info = paste0("incorrect dmvt calculation in R"))))

## test use through nimble function

nf <- nimbleFunction(
    run = function(x = double(1), mn = double(1),
                   scale = double(2), df = double(0)) {
        returnType(double(0))
        ch <- chol(scale)
        out <- dmvt_chol(x = x, mu = mn, cholesky = ch,
                         df = df, prec_param = FALSE, log = FALSE)
        return(out)
    }
)

cnf <- compileNimble(nf)

try(test_that("Test that dmvt_chol works correctly in R nimble function: ",
                  expect_equal(nf(x, mn, sc, df), (truth), 
                  info = paste0("incorrect dmvt value in R nimble function"))))

try(test_that("Test that dmvt_chol works correctly in compiled nimble function: ",
                  expect_equal(cnf(x, mn, sc, df), (truth), 
                  info = paste0("incorrect dmvt value in compiled nimble function"))))

## test use in model

mvt_code <- nimbleCode({
    x[1:3] ~ dmvt(mn[1:3], scale = sc[1:3,1:3], df = df)
})

mvt_model <- nimbleModel(mvt_code, constants = list(mn = mn, sc = sc, prec = FALSE, df = df))

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
# reference_samps <- mvtnorm::rmvt(n = 10000, delta = mn, sigma = sc, df = df)
r_samps <- t(replicate(10000, rmvt_chol(n = 1, mn, chol(sc), df, prec_param = FALSE)))
true_cov <- sc*df/(df-2)


try(test_that("Test that random samples (R) have correct mean: ",
              expect_equal(colMeans(r_samps), mn, 
                                                    tol = 0.03,
                          info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (R) have correct covariance: ",
              expect_equal(cov(r_samps), true_cov,
                                               tol = 0.1,
                          info = "Difference in covs exceeds tolerance")))


nf_sampling <- nimbleFunction(
    run = function(mn = double(1), scale = double(2), df = double(0)) {
        returnType(double(1))
        ch <- chol(scale)
        out <- rmvt_chol(n = 1, mu = mn, cholesky = ch,
                         df = df, prec_param = FALSE)
        return(out)
    }
)

nf_samps <- t(replicate(10000, nf_sampling(mn, sc, df)))

try(test_that("Test that random samples (nf) have correct mean: ",
              expect_equal(colMeans(nf_samps), mn, 
                                                    tol = 0.03,
                          info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (nf) have correct covariance: ",
              expect_equal(cov(nf_samps), true_cov, 
                                                tol = 0.1,
                          info = "Difference in covs exceeds tolerance")))

## sampling via `simulate`
simul_samp <- function(model) {
    model$simulate()
    return(model$x)
}

simul_samps <- t(replicate(10000, simul_samp(c_mvt_model)))

try(test_that("Test that random samples (simulate) have correct mean: ",
              expect_equal(colMeans(simul_samps), mn,
                                                    tol = 0.03,
                          info = "Difference in means exceeds tolerance")))

try(test_that("Test that random samples (simulate) have correct covariance: ",
              expect_equal(cov(simul_samps), true_cov, 
                                                tol = 0.1,
                          info = "Difference in covs exceeds tolerance")))


## dmulti

set.seed(0)
normGen <- rmulti(1, 1000, prob = c(.1, .1, .8))
set.seed(0)
unnormGen <- rmulti(1, 1000, prob = c(10, 10, 80))
normResult <- dmulti(normGen, prob = c(.1, .1, .8), log = TRUE)
unnormResult <- dmulti(normGen, prob = c(10, 10, 80), log = TRUE)

try(test_that("rmulti handles 'probs' that do not sum to one: ",
              expect_identical(normGen, unnormGen,
                          info = "normalized and unnormalized probabilities give different results")))
try(test_that("dmulti handles 'probs' that do not sum to one: ",
              expect_equal(normResult, unnormResult,
                          info = "normalized and unnormalized probabilities give different results")))
