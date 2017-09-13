source(system.file(file.path('tests', 'test_utils.R'), package = 'nimble'))

context("Testing of miscellaneous functionality.")

## Regression test for Issue #563.
test_that("while() works even when an intermediate variable is needed", {
    mynf = nimbleFunction(
    run=function(xi = double(1)) {
        returnType(double(0))
        newInd <- 1
        n <- length(xi)
        while(sum(xi == newInd) > 0 & newInd <= n){ 
            newInd <- newInd+1
        }
        return(newInd)
    })
    cmynf=compileNimble(mynf)
    expect_equal(mynf(c(1, 3, 4)),
                 cmynf(c(1, 3, 4)),
                 info = 'error in while() with intermediate-implemented expression')
})

test_that("copyExprClass works", {
    testExpr <- nimble:::RparseTree2ExprClasses(quote(a <- foo(bar(b, 3, c) + 2) + 4))
    testExprCopy <- nimble:::copyExprClass(testExpr)
    expect_identical(capture.output(testExpr$show()),
                     capture.output(testExprCopy$show()),
                     info = 'incorrect copying of exprClass object')
})

test_that("Test of full model check", {
    oldmsg <- geterrmessage()

    code <- nimbleCode({
        y ~ dnorm(mu, 1)
        mu ~ dnorm(0, 1)
    })
    try(m <- nimbleModel(code, data = list(y = 0), check = TRUE, name = 'test'))

    errmsg <- geterrmessage()    

    expect_equal(oldmsg, errmsg, info = "found error in running full model check")
})

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

test_that("literal double NaN", {
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

test_that("literal double NA", {
    nf <- nimbleFunction(
        run = function() {
            ans <- NA
            return(ans)
            returnType(double())
        })
    expect_true(is.na(nf()), 'NA uncompiled')

    # Known failure https://github.com/nimble-dev/nimble/issues/487
    expect_error({
        cnf <- compileNimble(nf)
        expect_true(is.na(cnf()), 'NA compiled')
    }, info = 'Literal NA is not supported in C++ code')
})

test_that('is.nan of scalar', {
    nf <- nimbleFunction(
        run = function(x = double(0)) {
            return(is.nan(x))
            returnType(logical())
        })
    cnf <- compileNimble(nf)
    for (x in c(0, 1, -Inf, Inf, NaN, -1.234)) {
        expect_identical(nf(x), cnf(x), info = paste('x =', x))
    }

    # Known failure https://github.com/nimble-dev/nimble/issues/487
    x <- NA
    expect_identical(nf(x), cnf(x), info = paste('x =', x))
})

test_that('is.na of scalar', {
    nf <- nimbleFunction(
        run = function(x = double(0)) {
            return(is.na(x))
            returnType(logical())
        })
    cnf <- compileNimble(nf)
    for (x in c(0, 1, -Inf, Inf, NA, -1.234)) {
        expect_identical(nf(x), cnf(x), info = paste('x =', x))
    }
    
    # Known failure https://github.com/nimble-dev/nimble/issues/487
    x <- NaN
    expect_failure(expect_identical(nf(x), cnf(x), info = paste('x =', x)))
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
