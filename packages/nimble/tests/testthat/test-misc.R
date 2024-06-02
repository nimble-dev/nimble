source(system.file(file.path('tests', 'testthat', 'test_utils.R'), package = 'nimble'))

RwarnLevel <- options('warn')$warn
options(warn = 1)
nimbleVerboseSetting <- nimbleOptions('verbose')
nimbleOptions(verbose = TRUE)

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
    expect_identical(capture.output(testExpr$print()),
                     capture.output(testExprCopy$print()),
                     info = 'incorrect copying of exprClass object')
})

test_that("Test of full model check", {
    oldmsg <- geterrmessage()

    code <- nimbleCode({
        y ~ dnorm(mu, 1)
        mu ~ dnorm(0, 1)
    })
    expect_message(m <- nimbleModel(code, data = list(y = 0), check = TRUE, name = 'test'),
                   "NAs were detected")

    errmsg <- geterrmessage()    

    expect_equal(oldmsg, errmsg, info = "found error in running full model check")
})

test_that("No nimKeyword appears in specificCallReplacements", {
    duplicates <- intersect(names(nimble:::nimKeyWords), names(nimble:::specificCallReplacements))
    expect_true(length(duplicates) == 0, "symbols appear in both nimKeywords and specificCallReplacements")
})

test_that("No nimKeyword appears in specificCallHandlers", {
    duplicates <- intersect(names(nimble:::nimKeyWords), names(nimble:::specificCallHandlers))
    expect_true(length(duplicates) == 0, "symbols appear in both nimKeywords and specificCallHandlers")
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

test_that("literal NA cases", {
    nf <- nimbleFunction(
        setup = TRUE,
        methods = list(
            ## doubles
            dvec = function(y = double(1)) {
                returnType(double(1))
                y[1] <- NA
                return(y)
            },
            dvec_add = function(y = double(1)) {
                returnType(double(1))
                y[1] <- y[1] + NA
                return(y)
            },
            dvec_mult = function(y = double(1)) {
                returnType(double(1))
                y[1] <- y[1] * NA
                return(y)
            },
            dscalar = function(y = double(0)) {
                returnType(double(0))
                y <- NA 
                return(y)
            },
            dscalar_input = function(y = double(0, default = as.numeric(NA))) {
                ## test will use y = NA
                returnType(logical(0))
                if(is.na(y)) ans <- TRUE else ans <- FALSE
                return(ans)
            },
            ## integers
            ## We are not sure that all cases of math with integer NA
            ## should work, but these basic tests work.
            ivec = function(y = integer(1)) {
                returnType(integer(1))
                y[1] <- NA
                return(y)
            },
            ivec_add = function(y = integer(1)) {
                returnType(integer(1))
                y[1] <- y[1] + NA
                return(y)
            },
            ivec_mult = function(y = integer(1)) {
                returnType(integer(1))
                y[1] <- y[1] * NA
                return(y)
            },
            iscalar = function(y = integer(0)) {
                returnType(integer(0))
                y <- NA 
                return(y)
            },
            iscalar_input = function(y = integer(0, default = as.integer(NA))) {
                returnType(logical(0))
                if(is.na(y)) ans <- TRUE else ans <- FALSE
                return(ans)
            },
            ## logical
            lvec_fails = function(y = logical(1)) {
                returnType(logical(1))
                y[1] <- NA ## I don't expect this to work
                return(y)
            },
            lscalar_fails = function(y = logical(0)) {
                returnType(logical(0))
                y <- NA ## I don't expect this to work
                return(y)
            },
            ## casting returns
            create_return_d = function() {
                ans <- NA
                return(ans)
                returnType(double())
            },
            create_return_i = function() {
                ans <- NA
                return(ans)
                returnType(integer())
            },
            create_return_l_fails = function() {
                ans <- NA
                return(ans)
                returnType(logical())
            }
        )
    )

    nf1 <- nf()
    cnf1 <- compileNimble(nf1)
    for(i in 1:2) {
        if(i == 1) {
            compiled <- FALSE
            f <- nf1
        } else {
            compiled <- TRUE
            f <- cnf1
        }
        expect_true(is.na(f$dvec(as.numeric(1:3))[1]), 'NA dvec')
        expect_true(is.na(f$dvec_add(as.numeric(1:3))[1]), 'NA dvec_add')
        expect_true(is.na(f$dvec_mult(as.numeric(1:3))[1]), 'NA dvec_mult')
        expect_true(is.na(f$dscalar(as.numeric(1))), 'NA dscalar')
        expect_true(f$dscalar_input(as.numeric(NA)), 'NA dscalar_input')
        ## KNOWN FAILURE:
        if(compiled) 
            expect_false(f$dscalar_input(NA), 'compiled NA dscalar_input fails to cast from logical')
        else
            expect_true(f$dscalar_input(NA), 'NA dscalar_input fails to cast from logical')
        
        expect_true(is.na(f$ivec(as.integer(1:3))[1]), 'NA ivec')
        expect_true(is.na(f$ivec_add(as.integer(1:3))[1]), 'NA ivec_add')
        expect_true(is.na(f$ivec_mult(as.integer(1:3))[1]), 'NA ivec_mult')
        expect_true(is.na(f$iscalar(as.integer(1))), 'NA iscalar')
        ## KNOWN FAILURE:
        if(compiled)
            expect_false(f$iscalar_input(as.integer(NA)), 'compiled NA dscalar_input fails')
        else
            expect_true(f$iscalar_input(as.integer(NA)), 'NA dscalar_input')
        ## KNOWN FAILURE:
        if(compiled)
            expect_false(f$iscalar_input(NA), 'compiled NA dscalar_input fails to cast from logical')
        else
            expect_true(f$iscalar_input(NA), 'NA dscalar_input')

        if(compiled)
            ## KNOWN FAILURE:
            expect_false(is.na(f$lscalar_fails(TRUE)), 'compiled NA lscalar_input fails')
        else
            expect_true(is.na(f$lscalar_fails(TRUE)), 'NA lscalar_input')

        if(compiled)
            ## KNOWN FAILURE:
            expect_false(is.na(f$lvec_fails(c(TRUE, TRUE))[1]), 'compiled NA lvec_input fails')
        else
            expect_true(is.na(f$lvec_fails(c(TRUE, TRUE))[1]), 'NA lvec_input')

        
        expect_true(is.na(f$create_return_d()), 'NA create_return_d')
        expect_true(is.na(f$create_return_i()), 'NA create_return_i')
        if(compiled)
            ## KNOWN FAILURE:
            expect_false(is.na(f$create_return_l_fails()), 'compiled NA create_return_l_fails fails')
        else
            expect_true(is.na(f$create_return_l_fails()), 'NA create_return_l_fails')
    }
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
    for (x in c(0, 1, -Inf, Inf, NA, NaN, -1.234)) {
        expect_identical(nf(x), cnf(x), info = paste('x =', x))
    }
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
})

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
})

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
})

## From GitHub Issue #695
test_that("setInits works in complicated case", {
    mc1 <- nimbleCode({
        for(i in 1:2) {
            for(j in f[i]:n) {
                x[i, j] ~ dnorm(0,1)
            }
        }
    })
    inits = list(x = matrix(1:6, nrow = 2))
    data = list(x = matrix(c(100, NA, 100, NA, NA, NA), nrow = 2))
    expect_message(m <- nimbleModel(mc1,
                     constants = list(f = c(1, 2),
                                      n = 3),
                     data = data,
                     inits = inits), "Ignoring non-NA values in `inits`")

    ## This model defines x[1, 1], x[1, 2], x[1, 3], x[2, 2], and x[2, 3]. 
    ## x[2,1] is not a node in the model.
    ## data are provided for x[1, 1] and x[1, 2]

    ##m$x ## This is wrong
    ##       [,1] [,2] [,3]
    ##[1,]   100   3    5
    ##[2,]   NA    4    NA
    ##matrix(c(100, 2, 100, 4, 5, 6), nrow = 2) ## It should be this
    ##       [,1] [,2] [,3]
    ##[1,]  100  100    5
    ##[2,]    2    4    6
    expect_equal(m$x, matrix(c(100, 2, 100, 4, 5, 6), nrow = 2))
    expect_identical(m$isDataEnv$'x', matrix(c(TRUE,FALSE,TRUE,FALSE,FALSE,FALSE),2))
    expect_length(m$isData('x[2,1]'), 0)
}
)

test_that("error trapping bad dimension size in setInits", {
    code <- nimbleCode({
        for(i in 1:3)
            for(j in 1:3)
                y[i,j]~dnorm(0,1)
    })
    m <- nimbleModel(code)
    expect_message(m$setInits(list(y = rnorm(3))), "Incorrect size or dimension of initial value")

    code <- nimbleCode({
        for(i in 1:3)
            for(j in 1:3)
                y[i,j]~dnorm(0,1)
    })
    m <- nimbleModel(code, data = list(y = matrix(rnorm(9),3)))
    expect_message(m$setInits(list(y = rnorm(3))), "Incorrect size or dimension of initial value")
    
    code <- nimbleCode({
        for(i in 1:3)
            for(j in 1:3)
                y[i,j]~dnorm(0,1)
    })
    expect_error(m <- nimbleModel(code, inits = list(y = rnorm(3))), "inconsistent dimensionality")
})

test_that("missing return statement is trapped",
{
    ## capture.output here just suppresses the print()ed error
    ## message to keep it out of Travis log files.
    expect_error(msg <- capture.output({a <- nimbleFunction(
                                            run = function() {
                                                returnType(double(1))
                                                out <- rnorm(10, 0, 1)
                                            })}),
                 info = "error-trapping missing return statement")
})

options(warn = RwarnLevel)
nimbleOptions(verbose = nimbleVerboseSetting)
