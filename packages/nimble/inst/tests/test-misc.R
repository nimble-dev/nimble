source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of miscellaneous functionality.")

oldmsg <- geterrmessage()

code <- nimbleCode({
    y ~ dnorm(mu, 1)
    mu ~ dnorm(0, 1)
})

try(m <- nimbleModel(code, data = list(y = 0), check = TRUE))

errmsg <- geterrmessage()    

test_that("Test of full model check", expect_equal(oldmsg, errmsg, info = "found error in running full model check"))

test_that("No nimKeyword appears in specificCallReplacements", {
    duplicates <- intersect(names(nimble:::nimKeyWords), names(nimble:::specificCallReplacements))
    if (length(duplicates) > 0) {
        fail(paste('These symbols appear in both nimKeywords and specificCallReplacements:', paste(duplicates, collapse = ', ')))
    }
})

test_that("No nimKeyword appears in specificCallHandlers", {
    duplicates <- intersect(names(nimble:::nimKeyWords), names(nimble:::specificCallHandlers))
    if (length(duplicates) > 0) {
        fail(paste('These symbols appear in both nimKeywords and specificCallHandlers:', paste(duplicates, collapse = ', ')))
    }
})

test_that("keyword next works", {
    nf <- nimbleFunction(
    run = function() {
        x <- numeric(2)
        for(i in 1:4) {
            if(i < 3) next
            x[i-2] <- i
        }
        return(x)
        returnType(double(1))
    },
    check = FALSE)
    nfRes <- nf()
    cnf <- compileNimble(nf)
    cnfRes <- cnf()
    expect_equal(nfRes, c(3, 4), info = 'keyword next uncompiled')
    expect_equal(cnfRes, c(3, 4), info = 'keyword next compiled')
})

test_that("literal NaN", {
    nf <- nimbleFunction(
        run = function() {
            ans <- NaN
            return(ans)
            returnType(double())
        })
    cnf <- compileNimble(nf)
    expect_true(is.nan(nf()), 'NaN uncompiled')
    expect_true(is.nan(cnf()), 'NaN compiled')
})

test_that("pi case 1", {
    nf <- nimbleFunction(
    run = function() {
        x = numeric(3)
        x[1] <- pi
        x[2] <- 2 * pi
        pi <- 10.1
        x[3] <- pi
        return(x)
        returnType(double(1))
    })
    cnf <- compileNimble(nf)
    expect_equal(nf(), c(pi, 2*pi, 10.1), info = 'pi case 1 uncompiled')
    expect_equal(cnf(), c(pi, 2*pi, 10.1), info = 'pi case 1 compiled')
    }
)

test_that("pi case 2", {
    nf <- nimbleFunction(
    run = function() {
        x = numeric(2)
        x[1] <- pi
        x[2] <- 2 * pi
        return(x)
        returnType(double(1))
    })
    cnf <- compileNimble(nf)
    expect_equal(nf(), c(pi, 2*pi), info = 'pi case 1 uncompiled')
    expect_equal(cnf(), c(pi, 2*pi), info = 'pi case 1 compiled')
    }
)

test_that("pi case 3", {
    nf <- nimbleFunction(
    run = function() {
        x = numeric(2)
        pi <- 10.1
        x[1] <- pi
        x[2] <- 2 * pi
        return(x)
        returnType(double(1))
    })
    cnf <- compileNimble(nf)
    expect_equal(nf(), c(10.1, 20.2), info = 'pi case 1 uncompiled')
    expect_equal(cnf(), c(10.1, 20.2), info = 'pi case 1 compiled')
    }
)
