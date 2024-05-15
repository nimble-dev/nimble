source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

test_that("Test that keyword processing for qbeta works", {
    nf <- nimbleFunction(
        run = function(p = double(0)) {
            returnType(double(0))
            out <- qbeta(shape2 = 1, shape1 = 1, log.p = TRUE, p = p)
            return(out)
        }
    )
    
    cnf <- compileNimble(nf)
    expect_equal(nf(log(0.5)), 0.5, 
                 info = paste0("incorrect processing of qbeta in R"))
    expect_equal(cnf(log(0.5)), 0.5, 
                 info = paste0("incorrect processing of qbeta in C"))
})


test_that("Test that keyword processing for rmnorm works", {
    nf <- nimbleFunction(
        run = function(mn = double(1), cov = double(2)) {
            returnType(double(1))
            ch <- chol(cov)
            out <- rmnorm_chol(cholesky = ch, prec_param = FALSE, mean = mn)
            return(out)
        }
    )
    cnf <- compileNimble(nf)
    
    mn <- c(1,2,3)
    cov <- diag(rep(1,3))
    cov[1,3] <- cov[3,1] <- 0.1
    
    set.seed(1)
    truth <- mn + t(chol(cov))%*%rnorm(3)
    attributes(truth) <- NULL  # otherwise compare matrix and vec

    set.seed(1)
    expect_equal(nf(mn, cov), truth, 
                 info = paste0("incorrect processing of rmnorm in R"))
    set.seed(1)
    expect_equal(cnf(mn, cov), truth, 
                 info = paste0("incorrect processing of rmnorm in C"))
})

# test pweibull as check of p to q names

test_that("Test that keyword processing for pweibull works", {
    nf <- nimbleFunction(
        run = function(q = double(0)) {
            returnType(double(0))
            out <- pweibull(lower.tail = 0, q = 0.5, scale = 1.5, shape = 3)
            return(out)
        }
    )
    
    cnf <- compileNimble(nf)
    
    q <- 0.5
    truth <- pweibull(q, 3, 1.5, lower.tail = FALSE)

    expect_equal(nf(q), truth, 
                 info = paste0("incorrect processing of pweibull in R"))
    expect_equal(cnf(q), truth, 
                 info = paste0("incorrect processing of pweibull in C"))
})

# dgamma involves special processing so check

test_that("Test that keyword processing for pgamma works", {
    nf <- nimbleFunction(
        run = function(q = double(0)) {
            returnType(double(0))
            out <- pgamma(lower.tail = 0, shape = 2, q = q, rate = 10)
            return(out)
        }
    )
    
    cnf <- compileNimble(nf)
    
    q <- 0.5
    truth <- pgamma(q, rate = 10, shape = 2, lower.tail = FALSE)
    expect_equal(nf(q), truth, 
                 info = paste0("incorrect processing of pgamma in R"))
    expect_equal(cnf(q), truth, 
                 info = paste0("incorrect processing of pgamma in C"))
})


test_that("Test that keyword processing for pgamma works", {
    nf <- nimbleFunction(
        run = function(q = double(0)) {
            returnType(double(0))
            out <- pgamma(lower.tail = 0, q = q, shape = 2, scale = 0.1)
            return(out)
        }
    )
    
    cnf <- compileNimble(nf)
    
    q <- 0.5
    truth <- pgamma(q, scale = 0.1, shape = 2, lower.tail = FALSE)
    expect_equal(nf(q), truth, 
                 info = paste0("incorrect processing of pgamma in R"))
    expect_equal(cnf(q), truth, 
                 info = paste0("incorrect processing of pgamma in C"))
})


# check both our dexp_nimble and base R dexp

test_that("Test that keyword processing for pexp works", {
    nf <- nimbleFunction(
        run = function(q = double(0)) {
            returnType(double(0))
            out <- pexp(lower.tail = 0, q = q, rate = 10)
            return(out)
        }
    )
    
    cnf <- compileNimble(nf)
    
    q <- 0.5
    truth <- pexp_nimble(q, rate = 10, lower.tail = FALSE)
    expect_equal(nf(q), truth, 
                 info = paste0("incorrect processing of pexp in R"))
    expect_equal(cnf(q), truth, 
                 info = paste0("incorrect processing of pexp in C"))
})

test_that("Test that keyword processing for pexp_nimble works", {
    nf <- nimbleFunction(
        run = function(q = double(0)) {
            returnType(double(0))
            out <- pexp_nimble(lower.tail = 0, q = q, rate = 10)
            return(out)
        }
    )
    
    cnf <- compileNimble(nf)
    
    q <- 0.5
    truth <- pexp_nimble(q, rate = 10, lower.tail = FALSE)
    expect_equal(nf(q), truth, 
                 info = paste0("incorrect processing of pexp_nimble in R"))
    expect_equal(cnf(q), truth, 
                 info = paste0("incorrect processing of pexp_nimble in C"))
})

test_that("Test that keyword processing for pexp_nimble with scale", {
    nf <- nimbleFunction(
        run = function(q = double(0)) {
            returnType(double(0))
            out <- pexp_nimble(lower.tail = 0, q = q, scale = 0.1)
            return(out)
        }
    )
    
    cnf <- compileNimble(nf)
    
    q <- 0.5
    truth <- pexp_nimble(q, scale = 0.1, lower.tail = FALSE)
    expect_equal(nf(q), truth, 
                 info = paste0("incorrect processing of pexp_nimble with scale in R"))
    expect_equal(cnf(q), truth, 
                 info = paste0("incorrect processing of pexp_nimble with scale in C"))
})
